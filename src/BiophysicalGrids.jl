module BiophysicalGrids

using Dates
using Statistics: mean

using Unitful

using Microclimate
using Microclimate: DEFAULT_DEPTHS, example_soil_moisture_model
using SolarRadiation
using FluidProperties
using FluidProperties: GoffGratch, Teten, Huang, VapourPressureEquation

using RasterDataSources
using NCDatasets       # triggers Rasters NCDatasets extension for NetCDF support
using Rasters
using Rasters: X, Y, Ti, Near

# ---------------------------------------------------------------------------
# Mesoclimate adjustments
# ---------------------------------------------------------------------------
include("Mesoclimate/Mesoclimate.jl")

# ---------------------------------------------------------------------------
# Weather data sources
# ---------------------------------------------------------------------------
include("WeatherDataSources/common.jl")
include("WeatherDataSources/TerraClimate.jl")

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
    # Cloud / radiation
    cloud_from_srad,
    # Wind
    wind_profile_adjust,
    # TerraClimate data source types (re-exported from RasterDataSources)
    TerraClimate,
    Historical,
    Plus2C,
    Plus4C,
    # Weather data
    get_weather,
    # Simulation
    simulate_microclimate,
    # FluidProperties vapour pressure methods (not exported by FluidProperties itself)
    VapourPressureEquation,
    GoffGratch,
    Teten,
    Huang

end
