# Example: microclimate simulation at a single point using ERA5 hourly
# reanalysis data, driven via `MicroVectorProblem`.

using MicroclimateMapper
using Microclimate
using RasterDataSources
using Microclimate: example_microclimate_problem, example_soil_profile
using Rasters
using Unitful
using Dates
using GeoInterface: Wrappers as GIW
using ZarrDatasets
#using Plots

ENV["CDSAPI_KEY"] = "your key"
ENV["RASTERDATASOURCES_PATH"] = "your folder"

points = [geocode("Mount Stirling, Victoria"),]#GIW.Point((135.0, -30.0))
dates = Date(2025, 1, 1):Day(1):Date(2025, 12, 31)

# `ERA5` (public, no CDS key) or `ECMWFERA5Land` (needs a CDS key, ~9km
# resolution where available -- see RasterDataSources._cds_credentials).
weather_source = ECMWFERA5Land

# ---------------------------------------------------------------------------
# Build the model
# ---------------------------------------------------------------------------

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
           20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0] .* u"cm"
heights = [0.01, 1.2]u"m"

micro_model = MicroModel(;
    depths,
    heights,
    soil_properties_model = CampbelldeVriesSoilProperties(;
        de_vries_shape_factor = 0.1,
        recirculation_power = 4.0,
        return_flow_threshold = 0.162,
    ),
    soil_hydraulic_model  = example_soil_hydraulic_model(),
    snow_model            = SnowModel(),
    config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
)

model = MicroMapModel(;
    micro_model,
    dem_source = SRTM,
    weather_source,
    init_source = ECMWFERA5Land, # uses ERA5 soil temperature and moisture as initial condition
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",    
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
# `using ZarrDatasets` is required at the top-level script to activate the
# Rasters Zarr backend for the ARCO-ERA5 store. Once the JSON version
# conflict between RasterDataSources/dev and ZarrDatasets is resolved this
# can become a direct dependency.
problem = MicroVectorProblem(; model, points, dates,
    soil_profile = example_soil_profile(depths),
)

@time output = solve(problem)

# plot(u"°C".(output.soil_temperature[point=1, depth=At(0.0)]))
# plot(output.snow_depth[point=1, Ti=Between(DateTime(2025, 7, 1), DateTime(2025, 8, 1))])
# plot(output.soil_moisture[point=1, depth=(0.05 .. 0.20)])
# plot(output.ground_dew[point=1, Ti=Between(DateTime(2025, 7, 1), DateTime(2025, 8, 1))])
# plot!(output.ground_frost[point=1, Ti=Between(DateTime(2025, 7, 1), DateTime(2025, 8, 1))])