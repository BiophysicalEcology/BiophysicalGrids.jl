# Gridded solar radiation across a mountainous terrain using a downloaded SRTM DEM.
#
# Study area: Mont Blanc massif, Chamonix valley, French Alps (~100 × 100 pixels)
# DEM source:  SRTM via RasterDataSources (CSI-CGIAR mirror, ~90 m / 3 arcsecond)
# Day:         Summer solstice (day 172) — maximises north/south-facing contrast
#
# Pipeline:
#   1. Download SRTM tile and crop to ~100 × 100 pixel study area
#   2. Reproject to UTM for metric slope/aspect/horizon calculations
#   3. Compute slope, aspect, and terrain horizon angles (Geomorphometry.jl)
#   4. Run SolarRadiation.jl per pixel for each hour
#   5. Integrate daily total and plot results (CairoMakie)
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
using Printf
using Plots
import Plots: heatmap, plot, savefig

# ============================================================================
# Configuration
# ============================================================================

# Study area: above Chamonix, where the Aiguilles Rouges face Mont Blanc.
# A ~100 × 100 pixel window at SRTM 3-arcsecond (~90 m) resolution.
# 100 pixels × 3" = 300" = 0.0833° lat  ≈ 9.3 km N–S
# 100 pixels × 3" / cos(45.9°) = 0.120° ≈ 8.6 km E–W
center_lon = 6.87    # °E
center_lat = 45.92   # °N
extent_lat = 0.0833  # degrees latitude
extent_lon = 0.120   # degrees longitude

lon_min = center_lon - extent_lon / 2
lon_max = center_lon + extent_lon / 2
lat_min = center_lat - extent_lat / 2
lat_max = center_lat + extent_lat / 2

simulation_day   = 172               # 21 June — summer solstice
hour_step        = 1.0
hours_of_day     = 6.0:hour_step:20.0
default_albedo   = 0.2
n_horizon_angles = 24

# ============================================================================
# Step 1: Download SRTM tile and crop to study area
# ============================================================================

println("Downloading SRTM tile covering study area...")
# getraster returns a matrix of tile paths (one per 5°×5° tile); the area
# fits within a single tile so we take the first valid path.
tile_paths  = getraster(SRTM; bounds = (lon_min, lat_min, lon_max, lat_max))
valid_paths = filter(!ismissing, vec(tile_paths))
isempty(valid_paths) && error("No SRTM tile found for the requested bounds.")
dem_full    = Raster(only(valid_paths))

println("Cropping to study area (~100 × 100 pixels)...")
dem_wgs84 = dem_full[X(Between(lon_min, lon_max)), Y(Between(lat_min, lat_max))]

nx_wgs = size(dem_wgs84, X)
ny_wgs = size(dem_wgs84, Y)
println("  Cropped size: $(nx_wgs) × $(ny_wgs) pixels")

# ============================================================================
# Step 2: Reproject to UTM for metric calculations
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
dem_data     = data_is_xy ? permutedims(dem_data_raw) : dem_data_raw  # → (Y, X)

# Geomorphometry needs (X, Y) with Y ascending (south → north).
dem_for_geomorph = data_is_xy ? dem_data_raw : permutedims(dem_data)
y_descending     = y_coords_utm[1] > y_coords_utm[end]
if y_descending
    dem_for_geomorph = reverse(dem_for_geomorph, dims = 2)
end

slope_g  = Geomorphometry.slope( dem_for_geomorph; method = Horn(), cellsize = cs)
aspect_g = Geomorphometry.aspect(dem_for_geomorph; method = Horn(), cellsize = cs)

slope_out  = y_descending ? reverse(slope_g,  dims = 2) : slope_g
aspect_out = y_descending ? reverse(aspect_g, dims = 2) : aspect_g

# Restore to raster's native dimension order
if data_is_xy
    utm_slope  = Raster(slope_out;  dims = dims(utm_dem))
    utm_aspect = Raster(aspect_out; dims = dims(utm_dem))
else
    utm_slope  = Raster(permutedims(slope_out);  dims = dims(utm_dem))
    utm_aspect = Raster(permutedims(aspect_out); dims = dims(utm_dem))
end

# ============================================================================
# Step 4: Build per-pixel lat/lon coordinate arrays
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
# Step 6: Prepare unit-tagged rasters
# ============================================================================

println("Preparing rasters with physical units...")

tag(r, u) = map(x -> (ismissing(x) || isnan(x)) ? missing : x * u, r)

elevation_m   = tag(utm_dem,    1.0u"m")
slope_deg     = tag(utm_slope,  1.0u"°")
aspect_deg    = tag(utm_aspect, 1.0u"°")
latitude_deg  = tag(lat_raster, 1.0u"°")
longitude_deg = tag(lon_raster, 1.0u"°")
albedo_r      = map(x -> ismissing(x) ? missing : default_albedo, elevation_m)
pressure_r    = map(e -> ismissing(e) ? missing : atmospheric_pressure(e), elevation_m)

solar_model = SolarProblem(; scattered_uv = false)

# ============================================================================
# Step 7: Compute solar radiation for each hour
# ============================================================================

println("Computing solar radiation — day $simulation_day " *
        "(hours $(first(hours_of_day))–$(last(hours_of_day)))...")

global_terrain_hours = Vector{Matrix}(undef, 0)

for hour in hours_of_day
    print("  hour $hour  \r")
    solar_grid = Matrix{Any}(undef, ny_utm, nx_utm)

    for j in 1:nx_utm, i in 1:ny_utm
        # Map loop indices → raster positional indices
        ri, rj = data_is_xy ? (j, i) : (i, j)

        lat  = latitude_deg[ri, rj];   lon  = longitude_deg[ri, rj]
        elev = elevation_m[ri,  rj];   slp  = slope_deg[ri,  rj]
        asp  = aspect_deg[ri,   rj];   alb  = albedo_r[ri,   rj]
        pres = pressure_r[ri,   rj]

        if any(ismissing.([lat, lon, elev, slp, asp, alb, pres]))
            solar_grid[i, j] = missing
            continue
        end

        terrain = SolarTerrain(;
            latitude             = lat,
            longitude            = lon,
            elevation            = elev,
            slope                = slp,
            aspect               = asp,
            albedo               = alb,
            atmospheric_pressure = pres,
            horizon_angles       = horizons_u[i, j, :],
        )

        solar_grid[i, j] = solar_radiation(
            solar_model;
            solar_terrain = terrain,
            days  = [Float64(simulation_day)],
            hours = [hour],
        )
    end

    push!(global_terrain_hours,
        map(c -> ismissing(c) ? missing : ustrip(c.global_terrain[1]), solar_grid))
end
println()

# ============================================================================
# Step 8: Daily-integrated radiation (trapezoidal rule)
# ============================================================================

println("Integrating daily total radiation...")
hours_vec   = collect(hours_of_day)
daily_Wh    = zeros(ny_utm, nx_utm)

for k in 1:(length(hours_vec) - 1)
    dt = hours_vec[k + 1] - hours_vec[k]
    for j in 1:nx_utm, i in 1:ny_utm
        v1, v2 = global_terrain_hours[k][i, j], global_terrain_hours[k + 1][i, j]
        (!ismissing(v1) && !ismissing(v2)) && (daily_Wh[i, j] += dt * (v1 + v2) / 2)
    end
end

daily_MJ = daily_Wh .* 0.0036  # Wh/m² → MJ/m²/day
valid_d  = filter(!iszero, vec(daily_MJ))
println("  Daily range: $(round(minimum(valid_d); digits=1)) – " *
        "$(round(maximum(valid_d); digits=1)) MJ/m²/day")

# ============================================================================
# Step 9: Plot results
# ============================================================================

println("Plotting...")
fmt_h(h) = @sprintf("%02d:%02d", floor(Int, h), round(Int, (h - floor(h)) * 60))

# Plots.heatmap requires y ascending (rasters are typically north-first / descending).
# slope_out / aspect_out come from Geomorphometry in (nx, ny) order; permute to (ny, nx).
y_plt = y_coords_utm[1] > y_coords_utm[end] ? reverse(y_coords_utm) : y_coords_utm
flipy(m) = y_coords_utm[1] > y_coords_utm[end] ? m[end:-1:1, :] : m

dem_plt    = flipy(dem_data)
slope_plt  = flipy(permutedims(slope_out))
aspect_plt = flipy(permutedims(aspect_out))

common_kw = (; aspect_ratio = :equal, xlabel = "Easting (m)", ylabel = "Northing (m)")

# ── Fig. 1: Terrain properties ─────────────────────────────────────────────
p_elev = heatmap(x_coords_utm, y_plt, dem_plt;
    color = :terrain, title = "Elevation (m)",
    clims = extrema(filter(!isnan, vec(dem_plt))), common_kw...)
p_slope = heatmap(x_coords_utm, y_plt, slope_plt;
    color = :YlOrRd, title = "Slope (°)",
    clims = (0.0, maximum(filter(!isnan, vec(slope_plt)))), common_kw...)
p_aspect = heatmap(x_coords_utm, y_plt, aspect_plt;
    color = :hsv, title = "Aspect (°)", clims = (0.0, 360.0), common_kw...)

display(plot(p_elev, p_slope, p_aspect;
    layout = (1, 3), size = (1200, 420), left_margin = 4Plots.mm,
    plot_title = "Terrain — Chamonix, French Alps  (SRTM ~90 m,  ~$(nx_utm)×$(ny_utm) pixels)"))
savefig("chamonix_terrain.png")
println("  Saved chamonix_terrain.png")

# ── Fig. 2: Horizon angles (4 cardinal directions) ─────────────────────────
hz_dirs  = [(1, "N 0°"), (7, "E 90°"), (13, "S 180°"), (19, "W 270°")]
hz_valid = filter(!isnan, vec(horizons))
hz_clims = isempty(hz_valid) ? (0.0, 1.0) : extrema(hz_valid)

hz_panels = [heatmap(x_coords_utm, y_plt, flipy(horizons[:, :, d]);
    color = :YlOrRd, title = "Horizon $lbl", clims = hz_clims,
    colorbar_title = "°", common_kw...) for (d, lbl) in hz_dirs]

display(plot(hz_panels...; layout = (1, 4), size = (1400, 420), left_margin = 4Plots.mm,
    plot_title = "Terrain horizon angles — Chamonix"))
savefig("chamonix_horizons.png")
println("  Saved chamonix_horizons.png")

# ── Fig. 3: Solar radiation — 2×2 hourly panels ────────────────────────────
all_vals  = vcat([Float64.(filter(!ismissing, vec(g))) for g in global_terrain_hours]...)
s_clims   = (0.0, maximum(all_vals))
panel_idx = round.(Int, range(2, length(hours_vec) - 1; length = 4))

solar_panels = [heatmap(x_coords_utm, y_plt,
    flipy(Float64.(coalesce.(global_terrain_hours[pi], NaN)));
    color = :inferno, title = "Hour $(fmt_h(hours_vec[pi]))",
    clims = s_clims, colorbar_title = "W/m²", common_kw...)
    for pi in panel_idx]

display(plot(solar_panels...; layout = (2, 2), size = (1100, 900), left_margin = 4Plots.mm,
    plot_title = "Solar radiation — Day $simulation_day (summer solstice), Chamonix"))
savefig("chamonix_solar_panel.png")
println("  Saved chamonix_solar_panel.png")

# ── Fig. 4: Daily total radiation ──────────────────────────────────────────
display(heatmap(x_coords_utm, y_plt, flipy(daily_MJ);
    color = :inferno, colorbar_title = "MJ/m²/day",
    title = "Daily total solar radiation — Day $simulation_day, Chamonix",
    size = (850, 720), left_margin = 6Plots.mm, common_kw...))
savefig("chamonix_solar_daily.png")
println("  Saved chamonix_solar_daily.png")

println("\nDone! $(nx_utm)×$(ny_utm) pixel grid, day $simulation_day, " *
        "$(length(hours_vec)) hours simulated.")
println("Daily solar range: $(round(minimum(valid_d); digits=1)) – " *
        "$(round(maximum(valid_d); digits=1)) MJ/m²/day")
