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
# Each pixel runs a single July representative day (24 hours) with up to iterate_day=10
# passes (the default), stopping early when the maximum nodal temperature change is
# below convergence_tolerance=0.1 K, so near-surface soil temperatures converge.
# Rows are processed sequentially; each pixel is warm-started from the converged
# midnight soil-temperature profile of the pixel directly above it, reducing the
# iterations needed to converge.  Columns within each row are parallelised.
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
using SolarRadiation
using FluidProperties
using FluidProperties: Teten, Huang, GoffGratch, VPLookupTable
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

# Vapour pressure formula used in humidity lapse correction and simulation.
# GoffGratch() is the most accurate but slowest; alternatives from FluidProperties:
#   Teten() — simple empirical, fastest
#   Huang()  — more accurate than Teten, faster than GoffGratch
#   VPLookupTable() — precomputed table, fastest for gridded runs
vp_method = VPLookupTable()


# Time snapshots: step index is 1-based (hour + 1)
snapshot_hours = collect(0:23)          # all 24 hours of the day
snapshot_steps = snapshot_hours .+ 1   # 1-based step index within the day
hour_labels    = [@sprintf("%02d:00", h) for h in snapshot_hours]
nhours         = length(snapshot_hours)   # 24

# Subset used for the static 2×3 panel plots (6 representative hours)
panel_ks  = [1, 7, 10, 13, 16, 19]    # indices into snapshot_hours: 00,06,09,12,15,18
panel_labels = ["Midnight", "Dawn", "Mid-morning", "Midday", "Mid-afternoon", "Dusk"]

# ============================================================================
# Steps 1–6: Load DEM, reproject to UTM, compute terrain grids
# ============================================================================

println("Downloading SRTM DEM and reprojecting to UTM...")
(; utm_dem, x_coords_utm, y_coords_utm, nx_utm, ny_utm, cs) =
    load_utm_dem(center_lon, center_lat, extent_lon, extent_lat)
println("  UTM grid: $(nx_utm) × $(ny_utm) pixels, " *
        "cell size ≈ $(round(cs[1]; digits=1)) × $(round(cs[2]; digits=1)) m")

println("Computing terrain grids (slope, aspect, horizons)...")
(; dem_data, data_is_xy, y_descending,
   elevation_m, slope_deg, aspect_deg,
   latitude_deg, longitude_deg, pressure_r,
   horizons_u) = compute_terrain_grids(utm_dem, x_coords_utm, y_coords_utm;
                                       n_horizon_angles)

# ============================================================================
# Step 7: TerraClimate weather — full year at center pixel, slice to July
#
# The full year is downloaded so that get_weather can compute the annual mean
# air temperature for use as the deep soil boundary condition. The July slice
# carries this annual mean in deep_soil_temperature (same value every month).
# ============================================================================

println("Obtaining TerraClimate weather for year $year...")

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

aerosol_optical_depth = get_aerosol_optical_depth(center_lat, center_lon, 0.01, july)
solar_model = SolarProblem(; aerosol_optical_depth)

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

# ---- Precompute per-pixel lapse-corrected weather and T_init values (single-threaded).
# Doing this outside the threaded loop keeps struct allocation and unit arithmetic
# out of the hot path, reducing GC pressure during the parallel simulation.
wp_grid     = Array{Any}(undef, ny_utm, nx_utm)   # lapse-corrected weather per pixel
Tmean_grid  = fill(NaN * u"K", ny_utm, nx_utm)    # July mean air temperature
Tdeep_grid  = fill(NaN * u"K", ny_utm, nx_utm)    # deep soil temperature

for I in CartesianIndices((ny_utm, nx_utm))
    i, j   = I[1], I[2]
    ri, rj = data_is_xy ? (j, i) : (i, j)
    elev   = elevation_m[ri, rj]
    ismissing(elev) && continue
    wp = lapse_correct_weather(weather_july, elev - center_elev_u)
    wp_grid[i, j]    = wp
    Tmean_grid[i, j] = (wp.environment_minmax.reference_temperature_min[1] +
                        wp.environment_minmax.reference_temperature_max[1]) / 2
    Tdeep_grid[i, j] = wp.environment_daily.deep_soil_temperature[1]
end

# Output arrays (ny, nx, nhours) — temperatures extracted in-loop; MicroResult discarded
T_soil0  = fill(NaN, ny_utm, nx_utm, nhours)  # soil surface (depth idx 1)
T_air1   = fill(NaN, ny_utm, nx_utm, nhours)  # air at 1 cm  (height idx 2)
T_air2   = fill(NaN, ny_utm, nx_utm, nhours)  # air at 2 m   (height idx 1, reference)
T_soil10 = fill(NaN, ny_utm, nx_utm, nhours)  # soil at 10 cm (depth idx 4)

n_total = nx_utm * ny_utm
n_done  = Threads.Atomic{Int}(0)

# Spatial warm-start: store the converged midnight soil-temperature profile for
# each completed pixel so the row below can use it as initial_soil_temperature.
# Rows are processed sequentially; columns within each row run in parallel.
# The first row uses an elevation-based estimate; subsequent rows use the
# converged profile from the pixel directly above (same column), which has
# similar elevation and hence a similar diurnal temperature regime.
T_converged = fill!(Matrix{Any}(undef, ny_utm, nx_utm), nothing)

@time for i in 1:ny_utm
    Threads.@threads :static for j in 1:nx_utm
        ri, rj = data_is_xy ? (j, i) : (i, j)

        elev = elevation_m[ri, rj]
        lat  = latitude_deg[ri, rj]
        lon  = longitude_deg[ri, rj]
        slp  = slope_deg[ri,    rj]
        asp  = aspect_deg[ri,   rj]
        pres = pressure_r[ri,   rj]

        if ismissing(elev) || ismissing(lat) || ismissing(lon) ||
           ismissing(slp)  || ismissing(asp) || ismissing(pres)
            continue
        end

        wp         = wp_grid[i, j]
        Tmean_july = Tmean_grid[i, j]
        Tdeep      = Tdeep_grid[i, j]

        # Warm-start from converged midnight profile of pixel above (if available).
        # Falls back to elevation-based estimate for the first row or missing neighbors.
        T_init = if i > 1 && !isnothing(T_converged[i - 1, j])
            T_converged[i - 1, j]
        else
            T_tmp      = Vector{typeof(1.0u"K")}(undef, 10)
            T_tmp[1:8] .= Tmean_july
            T_tmp[9]    = (Tmean_july + Tdeep) / 2
            T_tmp[10]   = Tdeep
            T_tmp
        end

        st = SolarTerrain(;
            latitude             = lat,
            longitude            = lon,
            elevation            = elev,
            slope                = slp,
            aspect               = asp,
            albedo               = 0.15,
            atmospheric_pressure = pres,
            horizon_angles       = @view(horizons_u[i, j, :]),
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
            vapour_pressure_equation = vp_method,
            iterate_day = 5,
        )

        # Save midnight (hour 0) profile as warm-start seed for the row below.
        T_converged[i, j] = collect(result.soil_temperature[1, :])

        @inbounds for (k, s) in enumerate(snapshot_steps)
            T_soil0[i,  j, k] = ustrip(u"°C", result.soil_temperature[s, 1])
            T_air1[i,   j, k] = ustrip(u"°C", result.profile[s].air_temperature[1])
            T_air2[i,   j, k] = ustrip(u"°C", result.profile[s].air_temperature[2])
            T_soil10[i, j, k] = ustrip(u"°C", result.soil_temperature[s, 4])
        end

        d = Threads.atomic_add!(n_done, 1)
        if Threads.threadid() == 1 && (d % nx_utm == 0 || d == n_total - 1)
            pct = round(100 * (d + 1) / n_total; digits = 1)
            print("  row $i/$ny_utm  ($(d+1)/$n_total, $pct%)   \r")
        end
    end
end
println("\nSimulation complete.")

# ============================================================================
# Step 10: Plots — one 2×3 figure per variable
# ============================================================================

y_plt = ascending_y(y_coords_utm, zeros(ny_utm, 1))[1]  # ascending y coords for Plots

common_kw = (; aspect_ratio = :equal, xlabel = "Easting (m)", ylabel = "Northing (m)")

function plot_variable(data4d, var_label, fname)
    all_vals = filter(!isnan, vec(data4d))
    isempty(all_vals) && return
    clims = (minimum(all_vals), maximum(all_vals))
    nframes = size(data4d, 3)

    # For a 2×3 layout pick 6 evenly-spaced frames; for ≤6 use all
    ks = nframes <= 6 ? (1:nframes) : panel_ks
    ls = nframes <= 6 ? hour_labels[1:nframes] : panel_labels

    panels = [heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, data4d[:, :, ks[n]])[2];
        color = cgrad(:RdYlBu, rev = true), clims = clims,
        title = ls[n], colorbar_title = "°C",
        titlefontsize = 9, common_kw...) for n in eachindex(ks)]

    display(plot(panels...; layout = (2, 3), size = (1400, 900),
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
        heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, data4d[:, :, k])[2];
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
plot_variable(T_air2,   "Air temperature at 2 m (°C)",    "chamonix_july_Tair2m.png")
plot_variable(T_air1,   "Air temperature at 1 cm (°C)",   "chamonix_july_Tair1cm.png")
plot_variable(T_soil0,  "Soil surface temperature (°C)",  "chamonix_july_Tsoil0.png")
plot_variable(T_soil10, "Soil temperature at 10 cm (°C)", "chamonix_july_Tsoil10cm.png")

println("Animating...")
animate_variable(T_air2,   "Air temperature at 2 m (°C)",    "chamonix_july_Tair2m.gif")
animate_variable(T_air1,   "Air temperature at 1 cm (°C)",   "chamonix_july_Tair1cm.gif")
animate_variable(T_soil0,  "Soil surface temperature (°C)",  "chamonix_july_Tsoil0.gif")
animate_variable(T_soil10, "Soil temperature at 10 cm (°C)", "chamonix_july_Tsoil10cm.gif")

println("\nDone. $(nx_utm)×$(ny_utm) pixel grid, July $year, " *
        "$(nhours) time snapshots ($(join(snapshot_hours, ", ")) h).")
