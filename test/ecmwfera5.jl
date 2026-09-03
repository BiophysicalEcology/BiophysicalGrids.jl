using Test
using Dates
using MicroclimateMapper
using MicroclimateMapper: variables, init_variables, native_group, canonical_name,
    fallback_source, fallback_layers, weather_calendar, native_timestep, loader,
    Daily, Hourly, ContiguousTimeSeries, MultiGroupContiguousTimeSeries,
    _contiguous_series_indices
using RasterDataSources
using Rasters.Extents: Extent

@testset "ECMWFERA5" begin
    @test weather_calendar(ECMWFERA5) == Daily()
    @test native_timestep(ECMWFERA5) == Hourly()
    @test loader(ECMWFERA5) == ContiguousTimeSeries()
    @test length(variables(ECMWFERA5)) == 9
    @test init_variables(ECMWFERA5) == variables(ECMWFERA5)  # no soil vars in :sfc
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

@testset "_contiguous_series_indices lon/lat conventions" begin
    # Lon 0-360 wraparound and descending latitude only hold for GCP
    # ARCO-ERA5, not ECMWF's own stores (ascending, -180..180). Cover both
    # with synthetic coords.
    epoch = DateTime(1970, 1, 1)
    hours = collect(0:23)  # one day, hourly
    time_start = epoch + Hour(0)
    time_end = epoch + Hour(5)

    # A point query is always buffered to a small range before it gets here
    # (see `_points_bbox`), never zero-width -- match that.
    buffered(lon, lat) = Extent(X = (lon - 0.5, lon + 0.5), Y = (lat - 0.5, lat + 0.5))

    # GCP ERA5: lon 0..360, lat descending 90..-90. Query Wisconsin
    # (-89.4557) -> must wrap to ~270.5 via ERA5's Longitude360 trait.
    gcp_coords = (; hours, epoch, lon = collect(0.0:1.0:359.0), lat = collect(90.0:-1.0:-90.0))
    r = _contiguous_series_indices(ERA5, gcp_coords, buffered(-89.4557, 43.0), time_start, time_end)
    @test all(x -> 270.0 <= x <= 271.0, r.xs)
    @test all(y -> 42.0 <= y <= 44.0, r.ys)
    @test length(r.ti) == 6

    # ECMWF stores: lon -180..180, lat ascending -90..90 -- default
    # (Longitude180) convention, no wraparound, no direction assumption.
    ecmwf_coords = (; hours, epoch, lon = collect(-180.0:1.0:179.0), lat = collect(-90.0:1.0:90.0))
    r2 = _contiguous_series_indices(ECMWFERA5Land, ecmwf_coords, buffered(-89.4557, 43.0), time_start, time_end)
    @test all(x -> -90.0 <= x <= -88.0, r2.xs)
    @test all(y -> 42.0 <= y <= 44.0, r2.ys)
    @test length(r2.ti) == 6
end
