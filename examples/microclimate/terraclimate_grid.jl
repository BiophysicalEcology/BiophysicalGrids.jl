# Gridded microclimate across the Chamonix SRTM DEM using TerraClimate forcing (July 2000).
#
# Pipeline:
#   1. Download SRTM DEM (~100×100 pixels) and reproject to UTM
#   2. Compute slope, aspect, horizon angles (Geomorphometry.jl)
#   3. Download TerraClimate weather (year 2000) at center pixel; slice to July
#   4. Per-pixel simulate_microclimate with lapse-rate-corrected weather
#   5. Plot soil surface T, air T at 1 cm, soil T at 5 cm, soil T at 20 cm
#      at 6 times of day: midnight, dawn, mid-morning, midday, mid-afternoon, dusk
#
# Each pixel runs a single July representative day (24 hours) with iterate_day=3
# (the default) so near-surface soil temperatures converge within the day.
# The deep soil boundary condition is the annual mean temperature at each pixel
# (lapse-corrected from center elevation), computed from the full TerraClimate
# annual record.
#
# Threading: start Julia with `julia --threads auto` for best performance.
#
# Extra packages beyond BiophysicalGrids core (install once):
#   using Pkg; Pkg.add(["ArchGDAL", "Plots"])

using BiophysicalGrids
using RasterDataSources
using Rasters, ArchGDAL
using Rasters.Lookups
using GeoFormatTypes: EPSG
using Geomorphometry
using SolarRadiation
using FluidProperties
using Unitful
using Statistics: median
using Printf
using Plots
import Plots: heatmap, plot, savefig

# ============================================================================
# Configuration
# ============================================================================

# Study area: above Chamonix, French Alps — same extent as grid_solar.jl
center_lon = 6.87     # °E
center_lat = 45.92    # °N
extent_lat = 0.0833   # ~100 SRTM pixels N–S
extent_lon = 0.120    # ~100 SRTM pixels E–W

lon_min = center_lon - extent_lon / 2
lon_max = center_lon + extent_lon / 2
lat_min = center_lat - extent_lat / 2
lat_max = center_lat + extent_lat / 2

year             = 2000
july             = 7          # month index
n_horizon_angles = 24

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 2.0]u"m"

# Vapour pressure formula used in humidity lapse correction.
# GoffGratch() is the most accurate but slowest; alternatives from FluidProperties:
#   Teten() — simple empirical, fastest
#   Huang()  — more accurate than Teten, faster than GoffGratch
vp_method = Teten()

# Time snapshots: step index is 1-based (hour + 1)
snapshot_hours = collect(0:23)          # all 24 hours of the day
snapshot_steps = snapshot_hours .+ 1   # 1-based step index within the day
hour_labels    = [@sprintf("%02d:00", h) for h in snapshot_hours]
nhours         = length(snapshot_hours)   # 24

# Subset used for the static 2×3 panel plots (6 representative hours)
panel_ks  = [1, 7, 10, 13, 16, 19]    # indices into snapshot_hours: 00,06,09,12,15,18
panel_labels = ["Midnight", "Dawn", "Mid-morning", "Midday", "Mid-afternoon", "Dusk"]

# ============================================================================
# Step 1: Download SRTM tile and crop to study area
# ============================================================================

println("Downloading SRTM tile covering study area...")
tile_paths  = getraster(SRTM; bounds = (lon_min, lat_min, lon_max, lat_max))
valid_paths = filter(!ismissing, vec(tile_paths))
isempty(valid_paths) && error("No SRTM tile found for the requested bounds.")
dem_full = Raster(only(valid_paths))

println("Cropping to study area (~100 × 100 pixels)...")
dem_wgs84 = dem_full[X(Between(lon_min, lon_max)), Y(Between(lat_min, lat_max))]
println("  Cropped size: $(size(dem_wgs84, X)) × $(size(dem_wgs84, Y)) pixels")

# ============================================================================
# Step 2: Reproject to UTM for metric slope/aspect/horizon calculations
# ============================================================================

println("Reprojecting to UTM...")
utm_crs = get_utm_crs(dem_wgs84)
utm_dem = Rasters.resample(dem_wgs84; crs = utm_crs, method = :bilinear)

x_coords_utm = collect(lookup(utm_dem, X))
y_coords_utm = collect(lookup(utm_dem, Y))
nx_utm       = length(x_coords_utm)
ny_utm       = length(y_coords_utm)
cs           = (abs(x_coords_utm[2] - x_coords_utm[1]),
                abs(y_coords_utm[2] - y_coords_utm[1]))

println("  UTM grid: $(nx_utm) × $(ny_utm) pixels, " *
        "cell size ≈ $(round(cs[1]; digits=1)) × $(round(cs[2]; digits=1)) m")

# ============================================================================
# Step 3: Slope and aspect via Geomorphometry.jl
# ============================================================================

println("Computing slope and aspect...")

# Rasters may store data in (X,Y) or (Y,X) internal order; detect which.
dem_data_raw = map(x -> ismissing(x) ? NaN : Float64(x), utm_dem.data)
data_is_xy   = size(dem_data_raw) == (nx_utm, ny_utm) && ny_utm != nx_utm
dem_data     = data_is_xy ? permutedims(dem_data_raw) : dem_data_raw  # → (ny, nx)

# Geomorphometry needs (nx, ny) with Y ascending (south → north).
dem_for_geomorph = data_is_xy ? dem_data_raw : permutedims(dem_data)
y_descending     = y_coords_utm[1] > y_coords_utm[end]
if y_descending
    dem_for_geomorph = reverse(dem_for_geomorph, dims = 2)
end

slope_g  = Geomorphometry.slope( dem_for_geomorph; method = Horn(), cellsize = cs)
aspect_g = Geomorphometry.aspect(dem_for_geomorph; method = Horn(), cellsize = cs)

slope_out  = y_descending ? reverse(slope_g,  dims = 2) : slope_g
aspect_out = y_descending ? reverse(aspect_g, dims = 2) : aspect_g

if data_is_xy
    utm_slope  = Raster(slope_out;  dims = dims(utm_dem))
    utm_aspect = Raster(aspect_out; dims = dims(utm_dem))
else
    utm_slope  = Raster(permutedims(slope_out);  dims = dims(utm_dem))
    utm_aspect = Raster(permutedims(aspect_out); dims = dims(utm_dem))
end

# ============================================================================
# Step 4: Per-pixel lat/lon coordinate arrays
# ============================================================================

println("Building lat/lon coordinate arrays...")

latlon_dem = Rasters.resample(utm_dem; crs = EPSG(4326))
lat_range  = lookup(latlon_dem, Y)
lon_range  = lookup(latlon_dem, X)
lat_south, lat_north = extrema(lat_range)
lon_west,  lon_east  = extrema(lon_range)

lats = y_coords_utm[1] < y_coords_utm[end] ?
    range(lat_south, lat_north; length = ny_utm) :
    range(lat_north, lat_south; length = ny_utm)
lons = x_coords_utm[1] < x_coords_utm[end] ?
    range(lon_west, lon_east; length = nx_utm) :
    range(lon_east, lon_west; length = nx_utm)

lat_matrix = repeat(collect(lats), 1, nx_utm)   # (ny, nx)
lon_matrix = repeat(collect(lons)', ny_utm, 1)  # (ny, nx)

if data_is_xy
    lat_raster = Raster(permutedims(lat_matrix), dims(utm_dem))
    lon_raster = Raster(permutedims(lon_matrix), dims(utm_dem))
else
    lat_raster = Raster(lat_matrix, dims(utm_dem))
    lon_raster = Raster(lon_matrix, dims(utm_dem))
end

# ============================================================================
# Step 5: Horizon angles from the DEM
# ============================================================================

println("Computing horizon angles ($n_horizon_angles directions)...")
horizons   = compute_horizon_angles(dem_data, x_coords_utm, y_coords_utm, n_horizon_angles)
horizons_u = horizons .* 1.0u"°"

# ============================================================================
# Step 6: Unit-tagged rasters
# ============================================================================

tag(r, u) = map(x -> (ismissing(x) || isnan(x)) ? missing : x * u, r)

elevation_m   = tag(utm_dem,    1.0u"m")
slope_deg     = tag(utm_slope,  1.0u"°")
aspect_deg    = tag(utm_aspect, 1.0u"°")
latitude_deg  = tag(lat_raster, 1.0u"°")
longitude_deg = tag(lon_raster, 1.0u"°")
pressure_r    = map(e -> ismissing(e) ? missing : atmospheric_pressure(e), elevation_m)

# ============================================================================
# Step 7: TerraClimate weather — full year at center pixel, slice to July
#
# The full year is downloaded so that get_weather can compute the annual mean
# air temperature for use as the deep soil boundary condition. The July slice
# carries this annual mean in deep_soil_temperature (same value every month).
# ============================================================================

println("Downloading TerraClimate weather for year $year...")

valid_elev    = filter(!isnan, vec(dem_data))
center_elev_u = median(valid_elev) * u"m"

weather = get_weather(TerraClimate, center_lon, center_lat;
    ystart    = year,
    elevation = center_elev_u,
)

function extract_month(ws, m)
    mm = ws.environment_minmax
    ed = ws.environment_daily
    eh = ws.environment_hourly
    new_mm = MonthlyMinMaxEnvironment(;
        reference_temperature_min = mm.reference_temperature_min[[m]],
        reference_temperature_max = mm.reference_temperature_max[[m]],
        reference_wind_min        = mm.reference_wind_min[[m]],
        reference_wind_max        = mm.reference_wind_max[[m]],
        reference_humidity_min    = mm.reference_humidity_min[[m]],
        reference_humidity_max    = mm.reference_humidity_max[[m]],
        cloud_min                 = mm.cloud_min[[m]],
        cloud_max                 = mm.cloud_max[[m]],
        minima_times              = mm.minima_times,
        maxima_times              = mm.maxima_times,
    )
    new_ed = DailyTimeseries(;
        shade                 = ed.shade[[m]],
        soil_wetness          = ed.soil_wetness[[m]],
        surface_emissivity    = ed.surface_emissivity[[m]],
        cloud_emissivity      = ed.cloud_emissivity[[m]],
        rainfall              = ed.rainfall[[m]],
        deep_soil_temperature = ed.deep_soil_temperature[[m]],
        leaf_area_index       = ed.leaf_area_index[[m]],
    )
    # Slice the hourly pressure vector: each month occupies 24 consecutive entries
    h_range = ((m - 1) * 24 + 1):(m * 24)
    new_eh = HourlyTimeseries(;
        pressure              = isnothing(eh.pressure)              ? nothing : eh.pressure[h_range],
        reference_temperature = isnothing(eh.reference_temperature) ? nothing : eh.reference_temperature[h_range],
        reference_humidity    = isnothing(eh.reference_humidity)    ? nothing : eh.reference_humidity[h_range],
        reference_wind_speed  = isnothing(eh.reference_wind_speed)  ? nothing : eh.reference_wind_speed[h_range],
        global_radiation      = isnothing(eh.global_radiation)      ? nothing : eh.global_radiation[h_range],
        longwave_radiation    = isnothing(eh.longwave_radiation)     ? nothing : eh.longwave_radiation[h_range],
        cloud_cover           = isnothing(eh.cloud_cover)           ? nothing : eh.cloud_cover[h_range],
        rainfall              = isnothing(eh.rainfall)              ? nothing : eh.rainfall[h_range],
        zenith_angle          = isnothing(eh.zenith_angle)          ? nothing : eh.zenith_angle[h_range],
    )
    return merge(ws, (;
        environment_minmax    = new_mm,
        environment_daily     = new_ed,
        environment_hourly    = new_eh,
        days                  = [ws.days[m]],
        soil_moisture_monthly = ws.soil_moisture_monthly[[m]],
    ))
end

weather_july = extract_month(weather, july)

center_elev_m = round(ustrip(u"m",  center_elev_u); digits = 0)
tmin_july_C   = round(ustrip(u"°C", weather_july.environment_minmax.reference_temperature_min[1]); digits = 1)
tmax_july_C   = round(ustrip(u"°C", weather_july.environment_minmax.reference_temperature_max[1]); digits = 1)
tdeep_C       = round(ustrip(u"°C", weather_july.environment_daily.deep_soil_temperature[1]);       digits = 1)
println("  Center elevation: $center_elev_m m")
println("  July Tmin: $tmin_july_C °C,  Tmax: $tmax_july_C °C,  deep soil T: $tdeep_C °C")

# ============================================================================
# Step 8: Shared soil model and lapse correction helper
# ============================================================================

soil_thermal_model = CampbelldeVriesSoilThermal(;
    bulk_density          = 1.3u"Mg/m^3",
    mineral_density       = 2.56u"Mg/m^3",
    de_vries_shape_factor = 0.1,
    mineral_conductivity  = 2.5u"W/m/K",
    mineral_heat_capacity = 870.0u"J/kg/K",
    saturation_moisture   = 0.42u"m^3/m^3",
    recirculation_power   = 4.0,
    return_flow_threshold = 0.162,
)

solar_model = SolarProblem()

# Lapse-correct temperature and humidity for a given elevation difference (Δz = pixel − center).
# Humidity is adjusted by conserving actual vapour pressure (dry adiabatic approximation)
# using BiophysicalGrids.rh_at_temperature.
function lapse_correct_weather(ws, elev_diff; method = vp_method)
    mm = ws.environment_minmax
    ed = ws.environment_daily
    T_min_new = lapse_adjust_temperature(mm.reference_temperature_min, elev_diff, EnvironmentalLapseRate())
    T_max_new = lapse_adjust_temperature(mm.reference_temperature_max, elev_diff, EnvironmentalLapseRate())
    new_mm = MonthlyMinMaxEnvironment(;
        reference_temperature_min = T_min_new,
        reference_temperature_max = T_max_new,
        reference_wind_min        = mm.reference_wind_min,
        reference_wind_max        = mm.reference_wind_max,
        # RH_min pairs with T_max time; RH_max pairs with T_min time
        reference_humidity_min    = rh_at_temperature(mm.reference_humidity_min, mm.reference_temperature_max, T_max_new, method),
        reference_humidity_max    = rh_at_temperature(mm.reference_humidity_max, mm.reference_temperature_min, T_min_new, method),
        cloud_min                 = mm.cloud_min,
        cloud_max                 = mm.cloud_max,
        minima_times              = mm.minima_times,
        maxima_times              = mm.maxima_times,
    )
    new_ed = DailyTimeseries(;
        shade                 = ed.shade,
        soil_wetness          = ed.soil_wetness,
        surface_emissivity    = ed.surface_emissivity,
        cloud_emissivity      = ed.cloud_emissivity,
        rainfall              = ed.rainfall,
        deep_soil_temperature = lapse_adjust_temperature(ed.deep_soil_temperature, elev_diff, EnvironmentalLapseRate()),
        leaf_area_index       = ed.leaf_area_index,
    )
    return merge(ws, (; environment_minmax = new_mm, environment_daily = new_ed))
end

# ============================================================================
# Step 9: Per-pixel microclimate simulation
# ============================================================================

println("Running per-pixel microclimate ($(nx_utm) × $(ny_utm) pixels, " *
        "$(Threads.nthreads()) thread(s))...")

# Output arrays: (ny, nx, nhours) — extract temperatures in-loop; discard MicroResult
T_soil0  = fill(NaN, ny_utm, nx_utm, nhours)  # soil surface (depth idx 1)
T_air1   = fill(NaN, ny_utm, nx_utm, nhours)  # air at 1 cm  (height idx 2, reversed)
T_soil5  = fill(NaN, ny_utm, nx_utm, nhours)  # soil at 5 cm (depth idx 3)
T_soil20 = fill(NaN, ny_utm, nx_utm, nhours)  # soil at 20 cm (depth idx 6)

progress_lock = ReentrantLock()
n_done        = Ref(0)
n_total       = nx_utm * ny_utm

Threads.@threads for j in 1:nx_utm
    for i in 1:ny_utm
        ri, rj = data_is_xy ? (j, i) : (i, j)

        elev = elevation_m[ri, rj]
        lat  = latitude_deg[ri,  rj]
        lon  = longitude_deg[ri, rj]
        slp  = slope_deg[ri,     rj]
        asp  = aspect_deg[ri,    rj]
        pres = pressure_r[ri,    rj]

        if any(ismissing.([elev, lat, lon, slp, asp, pres]))
            continue
        end

        elev_diff  = elev - center_elev_u
        wp         = lapse_correct_weather(weather_july, elev_diff)
        Tmean_july = (wp.environment_minmax.reference_temperature_min[1] +
                      wp.environment_minmax.reference_temperature_max[1]) / 2
        Tdeep      = wp.environment_daily.deep_soil_temperature[1]

        # Physically motivated initial soil temperature profile:
        # 0–50 cm (idx 1–8): July mean; 100 cm (idx 9): midpoint; 200 cm (idx 10): annual mean
        T_init = vcat(fill(Tmean_july, 8), [(Tmean_july + Tdeep) / 2], [Tdeep])

        st = SolarTerrain(;
            latitude             = lat,
            longitude            = lon,
            elevation            = elev,
            slope                = slp,
            aspect               = asp,
            albedo               = 0.15,
            atmospheric_pressure = pres,
            horizon_angles       = horizons_u[i, j, :],
        )
        mt = MicroTerrain(;
            elevation        = elev,
            roughness_height = 0.004u"m",
            karman_constant  = 0.4,
            dyer_constant    = 16.0,
            viewfactor       = 1.0,
        )

        result = simulate_microclimate(
            st, mt, soil_thermal_model, wp;
            depths, heights, solar_model,
            initial_soil_temperature = T_init,
        )

        for (k, s) in enumerate(snapshot_steps)
            T_soil0[i,  j, k] = ustrip(u"°C", result.soil_temperature[s, 1])
            T_air1[i,   j, k] = ustrip(u"°C", result.profile[s].air_temperature[2])
            T_soil5[i,  j, k] = ustrip(u"°C", result.soil_temperature[s, 3])
            T_soil20[i, j, k] = ustrip(u"°C", result.soil_temperature[s, 6])
        end

        lock(progress_lock) do
            n_done[] += 1
            if n_done[] % 200 == 0 || n_done[] == n_total
                pct = round(100 * n_done[] / n_total; digits = 1)
                print("  $(n_done[]) / $n_total ($pct%)   \r")
            end
        end
    end
end
println("\nSimulation complete.")

# ============================================================================
# Step 10: Plots — one 2×3 figure per variable
# ============================================================================

y_plt = y_descending ? reverse(y_coords_utm) : y_coords_utm
flipy(m) = y_descending ? m[end:-1:1, :] : m

common_kw = (; aspect_ratio = :equal, xlabel = "Easting (m)", ylabel = "Northing (m)")

function plot_variable(data4d, var_label, fname)
    all_vals = filter(!isnan, vec(data4d))
    isempty(all_vals) && return
    clims = (minimum(all_vals), maximum(all_vals))
    nframes = size(data4d, 3)

    # For a 2×3 layout pick 6 evenly-spaced frames; for ≤6 use all
    ks = nframes <= 6 ? (1:nframes) : panel_ks
    ls = nframes <= 6 ? hour_labels[1:nframes] : panel_labels

    panels = [heatmap(x_coords_utm, y_plt, flipy(data4d[:, :, ks[n]]);
        color = cgrad(:RdYlBu, rev = true), clims = clims,
        title = ls[n], colorbar_title = "°C",
        titlefontsize = 9, common_kw...) for n in eachindex(ks)]

    nr = nframes <= 6 ? 2 : 2
    nc = nframes <= 6 ? 3 : 3
    display(plot(panels...; layout = (nr, nc), size = (1400, 900),
        left_margin = 5Plots.mm,
        plot_title = "$var_label — Chamonix, July $year"))
    savefig(fname)
    println("  Saved $fname")
end

function animate_variable(data4d, var_label, fname; framerate = 4)
    all_vals = filter(!isnan, vec(data4d))
    isempty(all_vals) && return
    clims = (minimum(all_vals), maximum(all_vals))
    nframes = size(data4d, 3)
    # Labels: use snapshot_hours if lengths match, otherwise just frame indices
    labels_here = nframes == length(snapshot_hours) ?
        [@sprintf("%02d:00", snapshot_hours[k]) for k in 1:nframes] :
        [@sprintf("frame %d", k) for k in 1:nframes]

    anim = @animate for k in 1:nframes
        heatmap(x_coords_utm, y_plt, flipy(data4d[:, :, k]);
            color = cgrad(:RdYlBu, rev = true), clims = clims,
            title = "$var_label\n$(labels_here[k]) — Chamonix, July $year",
            xlabel = "Easting (m)", ylabel = "Northing (m)",
            colorbar_title = "°C", aspect_ratio = :equal,
            titlefontsize = 9, size = (700, 600),
            left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
    end
    gif(anim, fname; fps = framerate)
    println("  Saved $fname")
end

println("Plotting...")
plot_variable(T_soil0,  "Soil surface temperature (°C)",  "chamonix_july_Tsoil0.png")
plot_variable(T_air1,   "Air temperature at 1 cm (°C)",   "chamonix_july_Tair1cm.png")
plot_variable(T_soil5,  "Soil temperature at 5 cm (°C)",  "chamonix_july_Tsoil5cm.png")
plot_variable(T_soil20, "Soil temperature at 20 cm (°C)", "chamonix_july_Tsoil20cm.png")

println("Animating...")
animate_variable(T_soil0,  "Soil surface temperature (°C)",  "chamonix_july_Tsoil0.gif")
animate_variable(T_air1,   "Air temperature at 1 cm (°C)",   "chamonix_july_Tair1cm.gif")
animate_variable(T_soil5,  "Soil temperature at 5 cm (°C)",  "chamonix_july_Tsoil5cm.gif")
animate_variable(T_soil20, "Soil temperature at 20 cm (°C)", "chamonix_july_Tsoil20cm.gif")

println("\nDone. $(nx_utm)×$(ny_utm) pixel grid, July $year, " *
        "$(nhours) time snapshots ($(join(snapshot_hours, ", ")) h).")
