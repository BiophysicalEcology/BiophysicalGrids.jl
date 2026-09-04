using Test
using Dates
using MicroclimateMapper
using MicroclimateMapper: variables, init_variables, native_group, canonical_name,
    fallback_source, fallback_layers, weather_calendar, native_timestep, loader,
    Daily, Hourly, MultiGroupContiguousTimeSeries, _contiguous_series_indices
using RasterDataSources
using Rasters.Extents: Extent

@testset "ERA5ECMWF" begin
    @test weather_calendar(ERA5ECMWF) == Daily()
    @test native_timestep(ERA5ECMWF) == Hourly()
    @test loader(ERA5ECMWF) == MultiGroupContiguousTimeSeries()
    @test length(variables(ERA5ECMWF)) == 9
    @test init_variables(ERA5ECMWF) == variables(ERA5ECMWF)  # no soil vars in :sfc
    for var in variables(ERA5ECMWF)
        @test native_group(ERA5ECMWF, MicroclimateMapper.native_field(var)) == :sfc
    end
end

@testset "ERA5ECMWFLand" begin
    @test loader(ERA5ECMWFLand) == MultiGroupContiguousTimeSeries()
    @test fallback_source(ERA5ECMWFLand) == ERA5ECMWF
    @test fallback_layers(ERA5ECMWFLand) == (:cloud_cover,)

    # Every field in variables()/init_variables() must resolve to a group --
    # native_group errors (rather than silently returning something) for an
    # unknown field, so this only passes if the mapping is complete.
    for var in init_variables(ERA5ECMWFLand)
        field = MicroclimateMapper.native_field(var)
        @test native_group(ERA5ECMWFLand, field) isa Symbol
    end

    # Land + its ERA5ECMWF fallback together must cover every field ERA5
    # itself declares (the whole point of the fallback).
    land_names = Set(map(canonical_name, variables(ERA5ECMWFLand)))
    fallback_names = Set(fallback_layers(ERA5ECMWFLand))
    era5_names = Set(map(canonical_name, variables(ERA5ECMWF)))
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
    gcp_coords = (; hours, epoch, step_ms = 3_600_000, lon = collect(0.0:1.0:359.0), lat = collect(90.0:-1.0:-90.0))
    r = _contiguous_series_indices(ERA5, gcp_coords, buffered(-89.4557, 43.0), time_start, time_end)
    @test all(x -> 270.0 <= x <= 271.0, r.xs)
    @test all(y -> 42.0 <= y <= 44.0, r.ys)
    @test length(r.ti) == 6

    # ECMWF stores: lon -180..180, lat ascending -90..90 -- default
    # (Longitude180) convention, no wraparound, no direction assumption.
    ecmwf_coords = (; hours, epoch, step_ms = 3_600_000, lon = collect(-180.0:1.0:179.0), lat = collect(-90.0:1.0:90.0))
    r2 = _contiguous_series_indices(ERA5ECMWFLand, ecmwf_coords, buffered(-89.4557, 43.0), time_start, time_end)
    @test all(x -> -90.0 <= x <= -88.0, r2.xs)
    @test all(y -> 42.0 <= y <= 44.0, r2.ys)
    @test length(r2.ti) == 6

    # ERA5ECMWF's own :sfc group reports "seconds since ..." (confirmed
    # live), unlike ERA5-Land's "hours since ...". A hardcoded hour-step
    # would misread this by 3600x -- regression test for that bug.
    seconds = collect(0:3600:(23 * 3600))  # one day, hourly, but in seconds
    seconds_coords = (; hours = seconds, epoch, step_ms = 1_000,
        lon = collect(-180.0:1.0:179.0), lat = collect(-90.0:1.0:90.0))
    r3 = _contiguous_series_indices(ERA5ECMWF, seconds_coords, buffered(-89.4557, 43.0), time_start, time_end)
    @test length(r3.ti) == 6
end
