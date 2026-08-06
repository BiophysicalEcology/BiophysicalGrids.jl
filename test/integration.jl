# End-to-end raster solves over a tiny real DEM + weather grid, with and without
# lateral surface-water routing. Needs network access (SRTM + TerraClimate downloads),
# so these are integration tests rather than unit tests.
using Test
using MicroclimateMapper
using Microclimate: example_microclimate_problem, example_soil_profile, MicroConfig,
    DynamicSoilMoisture
using RasterDataSources, Rasters
using Rasters.Extents: Extent
using Dates, Unitful

# Routing needs a dynamic water balance (it infiltrates routed water); the example
# model defaults to PrescribedSoilMoisture, so build the inner model with
# DynamicSoilMoisture and reuse it for both runs.
micro_model = example_microclimate_problem(; config = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture())).model
soil_profile = example_soil_profile(micro_model.depths)
area = Extent(X = (146.00, 146.01), Y = (-36.00, -35.99))
dates = Date(2000, 1, 1):Day(1):Date(2000, 1, 2)

function _model(routing_model)
    MicroMapModel(; 
        micro_model,
        dem_source=SRTM,
        weather_source=TerraClimate{Historical},
        surface_albedo_source=0.15, 
        roughness_height_source=0.01u"m",
        routing_model
    )
end

base_layers = (:soil_temperature, :soil_moisture, :air_temperature,
    :relative_humidity, :wind_speed, :surface_water, :global_radiation,
    :sky_temperature, :snow_depth)

@testset "raster solve without routing" begin
    out = solve(MicroRasterProblem(; model=_model(nothing), area, dates, template=SRTM, soil_profile))
    for l in base_layers
        @test l in keys(out)
    end
    @test !(:runoff_generated in keys(out))     # no routing ⇒ no runoff layer
    @test ndims(out.soil_temperature) >= 3      # X × Y × Ti (× depth)
end

@testset "raster solve with routing" begin
    out = solve(MicroRasterProblem(; model=_model(SurfaceRunoffRouting()), area, dates, template=SRTM, soil_profile))
    for l in base_layers
        @test l in keys(out)                    # all base layers still produced
    end
    @test :runoff_generated in keys(out)        # routing adds the routed-runoff layer
    rg = out.runoff_generated
    @test size(rg)[1:2] == size(out.surface_water)[1:2]   # same spatial grid
    # Routed runoff is non-negative wherever it's finite (masked cells are NaN).
    @test all(x -> !isfinite(x) || x >= 0.0u"kg/m^2", skipmissing(rg))
end
