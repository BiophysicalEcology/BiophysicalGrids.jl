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

# ECMWFERA5/ECMWFERA5Land: authenticated CDS access, same coords/open shape
# but via RasterDataSources.open_zarr_array (injects the CDS bearer token
# Zarr.jl's own store can't send).
function MicroclimateMapper._contiguous_series_coords(source::RasterDataSources.CDSZarrSource)
    hours_arr = RasterDataSources.open_zarr_array(source, "time")
    lat_arr = RasterDataSources.open_zarr_array(source, "latitude")
    lon_arr = RasterDataSources.open_zarr_array(source, "longitude")
    (; hours = hours_arr[:], epoch = _parse_cf_epoch(hours_arr.attrs["units"]), lat = lat_arr[:], lon = lon_arr[:])
end
MicroclimateMapper._contiguous_series_open(source::RasterDataSources.CDSZarrSource, long_name::AbstractString) =
    RasterDataSources.open_zarr_array(source, long_name)

end
