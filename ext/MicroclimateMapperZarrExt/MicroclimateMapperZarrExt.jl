module MicroclimateMapperZarrExt

using MicroclimateMapper
using ZarrDatasets
using Dates
import RasterDataSources

const Zarr = ZarrDatasets.Zarr

_parse_cf_epoch(units) = DateTime(match(r"since (.+)$", units)[1], dateformat"yyyy-mm-dd HH:MM:SS")

function MicroclimateMapper._contiguous_series_coords(source::RasterDataSources.CachedCloudSource)
    hours_arr = Zarr.zopen(source.url * "/time")
    lat_arr = Zarr.zopen(source.url * "/latitude")
    lon_arr = Zarr.zopen(source.url * "/longitude")
    (; hours = hours_arr[:], epoch = _parse_cf_epoch(hours_arr.attrs["units"]), lat = lat_arr[:], lon = lon_arr[:])
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
    (; hours = hours_arr[:], epoch = _parse_cf_epoch(hours_arr.attrs["units"]), lat = lat_arr[:], lon = lon_arr[:])
end
MicroclimateMapper._contiguous_series_open(source::RasterDataSources.CDSZarrSource, long_name::AbstractString) =
    RasterDataSources.open_zarr_store(source).arrays[long_name]

end
