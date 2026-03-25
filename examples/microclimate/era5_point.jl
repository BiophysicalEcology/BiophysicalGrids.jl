# Example: microclimate simulation at a point using ERA5 hourly reanalysis data
#
# This script demonstrates the two-step API:
#   1. get_weather  — download ERA5 data and build environment objects
#   2. simulate_microclimate — assemble terrain/soil structs and run the model
#
# ERA5 provides actual hourly data for consecutive calendar days (not monthly
# typical days like TerraClimate). The returned DailyMinMaxEnvironment
# automatically enables daily=true mode in simulate_microclimate, so consecutive
# days inherit soil state and iterate once rather than converging independently.
#
# The location (-89.46°, 43.14°) matches the NicheMapR validation dataset so
# outputs can be compared against micro_era5 output.
#
# Extra packages (install once):
#   using Pkg; Pkg.add(["Plots", "CSV", "DataFrames"])

using BiophysicalGrids
using Microclimate
using SolarRadiation
using FluidProperties
using Unitful
using Dates

lon, lat  = -89.4557, 43.1379
elevation = 270.0u"m"

# ---------------------------------------------------------------------------
# Step 1: download and prepare ERA5 hourly weather forcing
# ---------------------------------------------------------------------------
# tstart/tend are inclusive DateTime bounds in UTC.
# grid_elevation: ERA5 model orography at the nearest 0.25° cell.  Omit for
# no lapse correction, or supply the ERA5 orography value for this cell
# (available from the Zarr store as the `z` variable / 9.80665 m).
println("Downloading ERA5 data...")
weather = get_weather(ERA5, lon, lat;
    tstart    = DateTime(2000, 1, 1),
    tend      = DateTime(2000, 12, 31, 23),
    elevation,
    # grid_elevation = 305.0u"m",         # ERA5 model orography if known
    vapour_pressure_method = GoffGratch(),
    lapse_rate_type        = EnvironmentalLapseRate(),
    use_era5_soil_moisture = false,        # set true to use ERA5 swvl1
)
println("  Got $(length(weather.days)) days of ERA5 data")

# ---------------------------------------------------------------------------
# Step 2: construct terrain structs (flat terrain for this example)
# ---------------------------------------------------------------------------
solar_terrain = SolarTerrain(;
    elevation,
    slope                = 0.0u"°",
    aspect               = 0.0u"°",
    horizon_angles       = fill(0.0u"°", 24),
    albedo               = 0.15,
    atmospheric_pressure = atmospheric_pressure(elevation),
    latitude             = lat * u"°",
    longitude            = lon * u"°",
)

micro_terrain = MicroTerrain(;
    elevation,
    roughness_height = 0.004u"m",
    karman_constant  = 0.4,
    dyer_constant    = 16.0,
    viewfactor       = 1.0,
)

# ---------------------------------------------------------------------------
# Step 3: construct soil thermal model
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Step 4: simulate
# ---------------------------------------------------------------------------
# daily=true is set automatically because weather.environment_minmax is a
# DailyMinMaxEnvironment — no need to specify it explicitly.
# Pass hourly_rainfall=true to use sub-daily (hourly) rainfall for soil
# moisture calculations instead of the daily totals.
println("Running microclimate simulation...")
result = simulate_microclimate(
    solar_terrain,
    micro_terrain,
    soil_thermal_model,
    weather;
    depths                   = [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]u"cm",
    heights                  = [0.01, 2.0]u"m",
    runmoist                 = false,
    clearsky                 = false,
    organic_soil_cap         = true,
    vapour_pressure_equation = GoffGratch(),
    # hourly_rainfall        = true,   # uncomment for sub-daily rainfall
)
println("  Done. Output size: $(size(result.soil_temperature))")  # (8784, ndepths) for year 2000

# Quick inspection of outputs
soil_T = result.soil_temperature         # Matrix (nhours × ndepths) in K
air_T  = [p.air_temperature for p in result.profile]  # per-height air temp per hour
@show size(soil_T)
@show result.soil_temperature[1, :]     # first hour, all depths

# ---------------------------------------------------------------------------
# Visual comparison (uncomment to run; requires Plots, CSV, DataFrames)
# ---------------------------------------------------------------------------
# using CSV, DataFrames, Plots
#
# let
#     data_dir = joinpath(@__DIR__, "..", "..", "test", "data", "micro_era5")
#     soil   = CSV.read(joinpath(data_dir, "soil_era5.csv"),   DataFrame)
#     metout = CSV.read(joinpath(data_dir, "metout_era5.csv"), DataFrame)
#
#     n = nrow(soil)
#     t = 1:n
#
#     # ---- NicheMapR reference vectors ----------------------------------------
#     soiltemps_nmr = Matrix(soil[:, ["D0cm","D2.5cm","D5cm","D10cm","D15cm","D20cm","D30cm","D50cm","D100cm","D200cm"]])
#     ta1cm_nmr  = (metout.TALOC .+ 273.15) .* 1u"K"
#     ta2m_nmr   = (metout.TAREF .+ 273.15) .* 1u"K"
#
#     # ---- Julia output -------------------------------------------------------
#     air_temperature_matrix = hcat([p.air_temperature for p in result.profile]...)'
#     depths_labels = ["$(round(ustrip(u"cm", d); digits=1)) cm"
#                      for d in [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]u"cm"]
#
#     p_st = plot(layout=(3, 3), size=(900, 800),
#                 title=reshape(depths_labels, 1, :), legend=:outertop)
#     for col in 1:9
#         plot!(p_st, t, ustrip.(u"°C", u"K".(result.soil_temperature[t, col]));
#               sp=col, label="Julia", color=:red, ylabel="°C")
#         plot!(p_st, t, soiltemps_nmr[:, col];
#               sp=col, label="NicheMapR", color=:black)
#     end
#     display(p_st)
#
#     p_air = plot(layout=(1, 2), size=(900, 400), legend=:outertop)
#     plot!(p_air, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 1]));
#           sp=1, label="Julia", color=:red, title="Air temp 1 cm", ylabel="°C")
#     plot!(p_air, t, ustrip.(u"°C", ta1cm_nmr);
#           sp=1, label="NicheMapR", color=:black)
#     plot!(p_air, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 2]));
#           sp=2, label="Julia", color=:red, title="Air temp 2 m")
#     plot!(p_air, t, ustrip.(u"°C", ta2m_nmr);
#           sp=2, label="NicheMapR", color=:black)
#     display(p_air)
# end
