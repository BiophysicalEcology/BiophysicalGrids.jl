using Test
using MicroclimateMapper
using MicroclimateMapper: variables, init_variables, native_group, canonical_name,
    fallback_source, fallback_layers, weather_calendar, native_timestep, loader,
    Daily, Hourly, ContiguousTimeSeries, MultiGroupContiguousTimeSeries
using RasterDataSources

@testset "ECMWFERA5" begin
    @test weather_calendar(ECMWFERA5) == Daily()
    @test native_timestep(ECMWFERA5) == Hourly()
    @test loader(ECMWFERA5) == ContiguousTimeSeries()
    @test length(variables(ECMWFERA5)) == 9
    @test length(init_variables(ECMWFERA5)) == 11  # + stl1, swvl1
end

@testset "ECMWFERA5Land" begin
    @test loader(ECMWFERA5Land) == MultiGroupContiguousTimeSeries()
    @test fallback_source(ECMWFERA5Land) == ECMWFERA5
    @test fallback_layers(ECMWFERA5Land) == (:cloud_cover,)

    # Every field in variables()/init_variables() must resolve to a group --
    # native_group errors (rather than silently returning something) for an
    # unknown field, so this only passes if the mapping is complete.
    for var in init_variables(ECMWFERA5Land)
        field = MicroclimateMapper.native_field(var)
        @test native_group(ECMWFERA5Land, field) isa Symbol
    end

    # Land + its ECMWFERA5 fallback together must cover every field ERA5
    # itself declares (the whole point of the fallback).
    land_names = Set(map(canonical_name, variables(ECMWFERA5Land)))
    fallback_names = Set(fallback_layers(ECMWFERA5Land))
    era5_names = Set(map(canonical_name, variables(ECMWFERA5)))
    @test era5_names ⊆ (land_names ∪ fallback_names)
end
