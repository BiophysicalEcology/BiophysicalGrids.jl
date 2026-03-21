# Example: microclimate simulation at a point using TerraClimate monthly data
#
# This script demonstrates the two-step API:
#   1. get_weather  — download TerraClimate data and build environment objects
#   2. simulate_microclimate — assemble terrain/soil structs and run the model
#
# The location (-89.46°, 43.14°) matches the NicheMapR validation dataset so
# outputs can be compared against micro_terra.R.
#
# For a real application, terrain parameters (elevation, slope, aspect,
# horizon angles) would typically be derived from a high-resolution DEM using
# Geomorphometry.jl:
#
#   using Geomorphometry, RasterDataSources
#   dem = getraster(SRTM30, ...)
#   elevation      = Geomorphometry.extract_elevation(dem, lon, lat)
#   slope, aspect  = Geomorphometry.slope_aspect(dem)
#   horizon_angles = Geomorphometry.horizon_angles(dem, lon, lat)
#
# Below we use flat terrain at a known elevation for simplicity.

using BiophysicalGrids
using Microclimate
using SolarRadiation
using FluidProperties
using Unitful

lon, lat = -89.4557, 43.1379
elevation = 270.0u"m"

# ---------------------------------------------------------------------------
# Step 1: download and prepare monthly weather forcing
# ---------------------------------------------------------------------------
# The lapse rate correction adjusts TerraClimate grid-cell temperatures
# to the site elevation. grid_elevation defaults to 0 m (sea level) when
# not supplied — provide it for better accuracy.
weather = get_weather(TerraClimate, lon, lat;
    ystart = 2000,
    elevation,
    # grid_elevation defaults to elevation (no lapse correction).
    # Provide the WorldClim 2.5-arcmin grid elevation at this cell for lapse correction,
    # e.g. grid_elevation = 308.0u"m" if WorldClim reports that for the cell.
    vapour_pressure_method = GoffGratch(),    # swap to Teten() or Huang()
    lapse_rate_type = EnvironmentalLapseRate(), # swap to DryAdiabaticLapseRate()
)

# ---------------------------------------------------------------------------
# Step 2: construct terrain structs
# ---------------------------------------------------------------------------
solar_terrain = SolarTerrain(;
    elevation,
    slope             = 0.0u"°",
    aspect            = 0.0u"°",
    horizon_angles    = fill(0.0u"°", 24),
    albedo            = 0.15,
    atmospheric_pressure = atmospheric_pressure(elevation),
    latitude          = lat * u"°",
    longitude         = lon * u"°",
)

micro_terrain = MicroTerrain(;
    elevation,
    roughness_height  = 0.004u"m",
    karman_constant   = 0.4,
    dyer_constant     = 16.0,
    viewfactor        = 1.0,
)

# ---------------------------------------------------------------------------
# Step 3: construct soil thermal model
# (CampbelldeVries parameters for a generic mineral soil)
# ---------------------------------------------------------------------------
soil_thermal_model = CampbelldeVriesSoilThermal(;
    bulk_density            = 1.3u"Mg/m^3",
    mineral_density         = 2.56u"Mg/m^3",
    de_vries_shape_factor   = 0.1,          # 0.33 for organic, 0.1 for mineral
    mineral_conductivity    = 2.5u"W/m/K",
    mineral_heat_capacity   = 870.0u"J/kg/K",
    saturation_moisture     = 0.42u"m^3/m^3",
    recirculation_power     = 4.0,
    return_flow_threshold   = 0.162,
)

# ---------------------------------------------------------------------------
# Step 4: optionally apply a climate change scenario
# ---------------------------------------------------------------------------
# Swap Historical for TerraClimate{Plus2C} or TerraClimate{Plus4C} to shift
# the weather forcing by the TerraClimate scenario deltas.  The baseline
# weather download (Step 1) is always from the historical record; the scenario
# only modifies the environment structs passed to simulate_microclimate.
weather_scenario = apply_climate_scenario(Historical, weather, lon, lat)
# weather_scenario = apply_climate_scenario(TerraClimate{Plus2C}, weather, lon, lat; ystart = 2000)
# weather_scenario = apply_climate_scenario(TerraClimate{Plus4C}, weather, lon, lat; ystart = 2000)

# ---------------------------------------------------------------------------
# Step 5: simulate
# ---------------------------------------------------------------------------
result = simulate_microclimate(
    solar_terrain,
    micro_terrain,
    soil_thermal_model,
    weather_scenario;
    depths   = [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]u"cm",
    heights  = [0.01, 2.0]u"m",
    runmoist = false,
    clearsky = false,
    organic_soil_cap = true,
)

# Quick inspection of outputs
soil_T = result.soil_temperature          # Matrix (timesteps × depths) in K
air_T  = [p.air_temperature for p in result.profile]  # per-height air temp
@show size(soil_T)
@show result.soil_temperature[1, :]      # first hour, all depths

# # ---------------------------------------------------------------------------
# # Sensitivity: swap vapour pressure equation
# # ---------------------------------------------------------------------------
# weather_teten = get_weather(TerraClimate, lon, lat;
#     ystart = 2000,
#     elevation,
#     vapour_pressure_method = Teten(),
# )
# result_teten = simulate_microclimate(
#     solar_terrain, micro_terrain, soil_thermal_model, weather_teten
# )

# # ---------------------------------------------------------------------------
# # Sensitivity: with dry adiabatic lapse rate
# # ---------------------------------------------------------------------------
# weather_dry = get_weather(TerraClimate, lon, lat;
#     ystart = 2000,
#     elevation,
#     lapse_rate_type = DryAdiabaticLapseRate(),
# )
# result_dry = simulate_microclimate(
#     solar_terrain, micro_terrain, soil_thermal_model, weather_dry
# )

# # ---------------------------------------------------------------------------
# # Multi-year run (2000–2002)
# # ---------------------------------------------------------------------------
# weather_3yr = get_weather(TerraClimate, lon, lat;
#     ystart = 2000, yfinish = 2002,
#     elevation,
# )
# result_3yr = simulate_microclimate(
#     solar_terrain, micro_terrain, soil_thermal_model, weather_3yr
# )
# @show size(result_3yr.soil_temperature)  # should be (36*24, ndepths)

# ---------------------------------------------------------------------------
# Visual comparisons against NicheMapR (run manually, not in CI)
# ---------------------------------------------------------------------------
# Requires Plots, CSV, DataFrames — install in your environment if needed:
#   using Pkg; Pkg.add(["Plots", "CSV", "DataFrames"])

using CSV, DataFrames, Plots

let
    data_dir = joinpath(@__DIR__, "..", "..", "test", "data", "micro_terra")
    soil   = CSV.read(joinpath(data_dir, "soil_monthly_terra.csv"),   DataFrame)
    metout = CSV.read(joinpath(data_dir, "metout_monthly_terra.csv"), DataFrame)

    n = nrow(soil)   # 288 = 12 months × 24 hours (year 2000)
    t = 1:n

    # ---- NicheMapR reference vectors ----------------------------------------
    soiltemps_nmr = Matrix(soil[:,   ["D0cm","D2.5cm","D5cm","D10cm","D15cm","D20cm","D30cm","D50cm","D100cm","D200cm"]])
    ta1cm_nmr  = (metout.TALOC .+ 273.15) .* 1u"K"
    ta2m_nmr   = (metout.TAREF .+ 273.15) .* 1u"K"
    rh1cm_nmr  =  metout.RHLOC ./ 100.0
    rh2m_nmr   =  metout.RH    ./ 100.0
    vel1cm_nmr =  metout.VLOC  .* 1u"m/s"
    vel2m_nmr  =  metout.VREF  .* 1u"m/s"

    # ---- Julia output matrices (ntimesteps × nheights) ----------------------
    air_temperature_matrix = hcat([p.air_temperature for p in result.profile]...)'
    humidity_matrix        = hcat([p.relative_humidity for p in result.profile]...)'
    wind_matrix            = hcat([p.wind_speed for p in result.profile]...)'

    depths_labels = ["$(round(ustrip(u"cm", d); digits=1)) cm"
                     for d in [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]u"cm"]

    # ---- Soil temperatures --------------------------------------------------
    soil0_julia = ustrip.(u"°C", u"K".(result.soil_temperature[t, 1]))
    soil0_nmr   = soiltemps_nmr[:, 1]
    soil_ylim   = (floor(min(minimum(soil0_julia), minimum(soil0_nmr))) - 1,
                   ceil( max(maximum(soil0_julia), maximum(soil0_nmr))) + 1)

    p_st = plot(layout=(3, 3), size=(900, 800),
                title=reshape(depths_labels, 1, :), legend=:outertop)
    for col in 1:9
        plot!(p_st, t, ustrip.(u"°C", u"K".(result.soil_temperature[t, col]));
              sp=col, label="Julia",     color=:red,   ylabel="°C", ylims=soil_ylim)
        plot!(p_st, t, soiltemps_nmr[:, col];
              sp=col, label="NicheMapR", color=:black)
    end
    display(p_st)

    # ---- Atmospheric profiles -----------------------------------------------
    p_atm = plot(layout=(3, 2), size=(900, 800), legend=:outertop)

    plot!(p_atm, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 1]));
          sp=1, label="Julia", color=:red, title="Air temp 1 cm", ylabel="°C")
    plot!(p_atm, t, ustrip.(u"°C", ta1cm_nmr);
          sp=1, label="NicheMapR", color=:black)

    plot!(p_atm, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 2]));
          sp=2, label="Julia", color=:red, title="Air temp 2 m")
    plot!(p_atm, t, ustrip.(u"°C", ta2m_nmr);
          sp=2, label="NicheMapR", color=:black)

    plot!(p_atm, t, humidity_matrix[t, 1];
          sp=3, label="Julia", color=:red, title="RH 1 cm", ylabel="–")
    plot!(p_atm, t, rh1cm_nmr;
          sp=3, label="NicheMapR", color=:black)

    plot!(p_atm, t, humidity_matrix[t, 2];
          sp=4, label="Julia", color=:red, title="RH 2 m")
    plot!(p_atm, t, rh2m_nmr;
          sp=4, label="NicheMapR", color=:black)

    plot!(p_atm, t, ustrip.(u"m/s", wind_matrix[t, 1]);
          sp=5, label="Julia", color=:red, title="Wind 1 cm", ylabel="m/s")
    plot!(p_atm, t, ustrip.(u"m/s", vel1cm_nmr);
          sp=5, label="NicheMapR", color=:black)

    plot!(p_atm, t, ustrip.(u"m/s", wind_matrix[t, 2]);
          sp=6, label="Julia", color=:red, title="Wind 2 m")
    plot!(p_atm, t, ustrip.(u"m/s", vel2m_nmr);
          sp=6, label="NicheMapR", color=:black)

    display(p_atm)
end
