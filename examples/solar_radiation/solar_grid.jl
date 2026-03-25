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
n_horizon_angles = 32

# ============================================================================
# Step 1: Download SRTM tile and crop to study area
# ============================================================================

println("Downloading SRTM DEM and reprojecting to UTM...")
(; utm_dem, x_coords_utm, y_coords_utm, nx_utm, ny_utm, cs) =
    load_utm_dem(center_lon, center_lat, extent_lon, extent_lat)
println("  UTM grid: $(nx_utm) × $(ny_utm) pixels, " *
        "cell size ≈ $(round(cs[1]; digits=1)) × $(round(cs[2]; digits=1)) m")

# ============================================================================
# Steps 2–5: Slope, aspect, lat/lon, horizon angles, unit-tagged rasters
# ============================================================================

println("Computing terrain grids (slope, aspect, horizons)...")
(; dem_data, data_is_xy, y_descending,
   elevation_m, slope_deg, aspect_deg,
   latitude_deg, longitude_deg, pressure_r,
   horizons_u) = compute_terrain_grids(utm_dem, x_coords_utm, y_coords_utm;
                                       n_horizon_angles)

albedo_r = map(x -> ismissing(x) ? missing : default_albedo, elevation_m)

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

# ascending_y(y, m) flips both y and the (ny,nx) matrix so heatmap gets ascending y.
# Slope/aspect rasters are unit-tagged; strip units and reorder to (ny, nx) first.
to_mat(r) = map(x -> ismissing(x) ? NaN : ustrip(x), r)
to_ny_nx(r) = data_is_xy ? permutedims(to_mat(r)) : to_mat(r)

y_plt, dem_plt    = ascending_y(y_coords_utm, dem_data)
_,     slope_plt  = ascending_y(y_coords_utm, to_ny_nx(slope_deg))
_,     aspect_plt = ascending_y(y_coords_utm, to_ny_nx(aspect_deg))

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
horizons = ustrip.(u"°", horizons_u)  # raw Float64 for plotting
hz_valid = filter(!isnan, vec(horizons))
hz_clims = isempty(hz_valid) ? (0.0, 1.0) : extrema(hz_valid)

hz_panels = [heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, horizons[:, :, d])[2];
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
    ascending_y(y_coords_utm, Float64.(coalesce.(global_terrain_hours[pi], NaN)))[2];
    color = :inferno, title = "Hour $(fmt_h(hours_vec[pi]))",
    clims = s_clims, colorbar_title = "W/m²", common_kw...)
    for pi in panel_idx]

display(plot(solar_panels...; layout = (2, 2), size = (1100, 900), left_margin = 4Plots.mm,
    plot_title = "Solar radiation — Day $simulation_day (summer solstice), Chamonix"))
savefig("chamonix_solar_panel.png")
println("  Saved chamonix_solar_panel.png")

# ── Fig. 4: Daily total radiation ──────────────────────────────────────────
display(heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, daily_MJ)[2];
    color = :inferno, colorbar_title = "MJ/m²/day",
    title = "Daily total solar radiation — Day $simulation_day, Chamonix",
    size = (850, 720), left_margin = 6Plots.mm, common_kw...))
savefig("chamonix_solar_daily.png")
println("  Saved chamonix_solar_daily.png")

# ── Fig. 5: Animated solar radiation — one frame per hour ──────────────────
println("Animating...")

solar_clims = (0.0, maximum(all_vals))

anim_solar = @animate for k in 1:length(hours_vec)
    frame = Float64.(coalesce.(global_terrain_hours[k], NaN))
    _, frame_plt = ascending_y(y_coords_utm, frame)
    heatmap(x_coords_utm, y_plt, frame_plt;
        color = :inferno, clims = solar_clims,
        title = "Global terrain radiation — $(fmt_h(hours_vec[k]))\nDay $simulation_day (summer solstice), Chamonix",
        colorbar_title = "W/m²", common_kw...,
        titlefontsize = 9, size = (700, 600),
        left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
end
gif(anim_solar, "chamonix_solar.gif"; fps = 4)
println("  Saved chamonix_solar.gif")

println("\nDone! $(nx_utm)×$(ny_utm) pixel grid, day $simulation_day, " *
        "$(length(hours_vec)) hours simulated.")
println("Daily solar range: $(round(minimum(valid_d); digits=1)) – " *
        "$(round(maximum(valid_d); digits=1)) MJ/m²/day")
