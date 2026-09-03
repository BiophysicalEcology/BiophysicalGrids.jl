module MicroclimateMapperZarrExt

using MicroclimateMapper
using ZarrDatasets
using Dates
import RasterDataSources

const Zarr = ZarrDatasets.Zarr

# Step unit varies by store -- ERA5-Land uses "hours since ...", ECMWFERA5's
# own :sfc group uses "seconds since ..." (confirmed live); hardcoding hours
# silently misreads one as the other.
const _CF_TIME_STEP_MS = Dict("hours" => 3_600_000, "days" => 86_400_000, "minutes" => 60_000, "seconds" => 1_000)
function _parse_cf_time(units)
    step_name, epoch_str = match(r"^(hours|days|seconds|minutes) since (.+)$", units).captures
    (; epoch = DateTime(epoch_str, dateformat"yyyy-mm-dd HH:MM:SS"), step_ms = _CF_TIME_STEP_MS[step_name])
end

function MicroclimateMapper._contiguous_series_coords(source::RasterDataSources.CachedCloudSource)
    hours_arr = Zarr.zopen(source.url * "/time")
    lat_arr = Zarr.zopen(source.url * "/latitude")
    lon_arr = Zarr.zopen(source.url * "/longitude")
    (; hours = hours_arr[:], _parse_cf_time(hours_arr.attrs["units"])..., lat = lat_arr[:], lon = lon_arr[:])
end

# `RasterStack(url; source=Zarrsource())`'s consolidated-metadata discovery
# only surfaces a subset of this store's ~280 arrays, so named variables are
# opened directly by their store subpath instead.
MicroclimateMapper._contiguous_series_open(source::RasterDataSources.CachedCloudSource, long_name::AbstractString) =
    Zarr.zopen(source.url * "/" * long_name)

# ECMWFERA5/ECMWFERA5Land: authenticated CDS access via
# RasterDataSources.open_zarr_store, which disk-caches chunks (open_zarr_array
# doesn't -- fine for a one-off array but far too slow for a full time series).
function MicroclimateMapper._contiguous_series_coords(source::RasterDataSources.CDSZarrSource)
    ds = RasterDataSources.open_zarr_store(source)
    hours_arr, lat_arr, lon_arr = ds.arrays["time"], ds.arrays["latitude"], ds.arrays["longitude"]
    (; hours = hours_arr[:], _parse_cf_time(hours_arr.attrs["units"])..., lat = lat_arr[:], lon = lon_arr[:])
end
MicroclimateMapper._contiguous_series_open(source::RasterDataSources.CDSZarrSource, long_name::AbstractString) =
    RasterDataSources.open_zarr_store(source).arrays[long_name]

end
