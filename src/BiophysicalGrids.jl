module BiophysicalGrids

using Dates
using Statistics: mean

using Unitful

using Microclimate
using Microclimate: DEFAULT_DEPTHS, example_soil_moisture_model
using SolarRadiation
using FluidProperties
using FluidProperties: GoffGratch, Teten, Huang, VapourPressureEquation, VPLookupTable

using GeoFormatTypes
using GeoInterface
using RasterDataSources
using RasterDataSources: layername
using DimensionalData: DimArray, Dim, At, dims
using NCDatasets       # triggers Rasters NCDatasets extension for NetCDF support
using Zarr
using Rasters
using Rasters: X, Y, Ti, Near, Between, lookup
using Geomorphometry
using Geomorphometry: Horn

const Point = GeoInterface.Wrappers.Point

"""longitude(point) → degrees east — extract longitude from a GeoInterface point."""
longitude(point) = GeoInterface.x(point)

"""latitude(point) → degrees north — extract latitude from a GeoInterface point."""
latitude(point) = GeoInterface.y(point)

# ---------------------------------------------------------------------------
# Atmosphere
# ---------------------------------------------------------------------------
include("Atmosphere/aerosol.jl")

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
# Grid simulation
# ---------------------------------------------------------------------------
include("Grid/solar_grid.jl")
include("Grid/microclimate_grid.jl")

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

export
    # Point / extent geometry
    Point,
    longitude,
    latitude,
    Extent,
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
    # Grid simulation
    solar_radiation_grid,
    GridInitStrategy,
    RowWarmStart,
    ERA5SoilInit,
    SolarMatchInit,
    LatitudeMatchInit,
    simulate_microclimate_grid,
    # Aerosol
    get_aerosol_optical_depth,
    # Terrain utilities
    get_utm_crs,
    load_utm_dem,
    compute_terrain_grids,
    compute_horizon_angles,
    ascending_y,
    # FluidProperties vapour pressure methods
    VapourPressureEquation,
    GoffGratch,
    Teten,
    Huang,
    VPLookupTable

end
