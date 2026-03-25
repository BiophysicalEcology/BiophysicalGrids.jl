module BiophysicalGrids

using Dates
using Statistics: mean

using Unitful

using Microclimate
using Microclimate: DEFAULT_DEPTHS, example_soil_moisture_model
using SolarRadiation
using FluidProperties
using FluidProperties: GoffGratch, Teten, Huang, VapourPressureEquation

using GeoFormatTypes
using RasterDataSources
using RasterDataSources: layername
using NCDatasets       # triggers Rasters NCDatasets extension for NetCDF support
using Zarr
using Rasters
using Rasters: X, Y, Ti, Near, Between, lookup
using Geomorphometry
using Geomorphometry: Horn

# ---------------------------------------------------------------------------
# Terrain utilities
# ---------------------------------------------------------------------------
include("Terrain/terrain_utils.jl")

# ---------------------------------------------------------------------------
# Mesoclimate adjustments
# ---------------------------------------------------------------------------
include("Mesoclimate/Mesoclimate.jl")

# ---------------------------------------------------------------------------
# Weather data sources
# ---------------------------------------------------------------------------
include("WeatherDataSources/common.jl")
include("WeatherDataSources/TerraClimate.jl")
include("WeatherDataSources/ERA5.jl")
include("WeatherDataSources/climate_scenarios.jl")

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

export
    # Lapse rate types
    LapseRate,
    EnvironmentalLapseRate,
    DryAdiabaticLapseRate,
    SaturatedAdiabaticLapseRate,
    CustomLapseRate,
    # Lapse rate functions
    lapse_rate,
    lapse_adjust_temperature,
    dewpoint_lapse_adjust,
    rh_at_temperature,
    # Cloud / radiation
    cloud_from_srad,
    # Wind
    wind_profile_adjust,
    # Data source types (re-exported from RasterDataSources)
    TerraClimate,
    Historical,
    Plus2C,
    Plus4C,
    ERA5,
    # Weather data
    get_weather,
    apply_climate_scenario,
    # Microclimate types (re-exported from Microclimate.jl)
    MicroTerrain,
    CampbelldeVriesSoilThermal,
    DailyTimeseries,
    MonthlyMinMaxEnvironment,
    DailyMinMaxEnvironment,
    HourlyTimeseries,
    # Simulation
    simulate_microclimate,
    # Terrain utilities
    get_utm_crs,
    load_utm_dem,
    compute_terrain_grids,
    compute_horizon_angles,
    ascending_y,
    # FluidProperties vapour pressure methods (not exported by FluidProperties itself)
    VapourPressureEquation,
    GoffGratch,
    Teten,
    Huang

end
