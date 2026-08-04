module MicroclimateMapperZarrExt

using MicroclimateMapper
using ZarrDatasets
using Dates

const Zarr = ZarrDatasets.Zarr

function MicroclimateMapper._contiguous_series_coords(url::AbstractString)
    hours_arr = Zarr.zopen(url * "/time")
    lat_arr = Zarr.zopen(url * "/latitude")
    lon_arr = Zarr.zopen(url * "/longitude")
    epoch = DateTime(match(r"since (.+)$", hours_arr.attrs["units"])[1], dateformat"yyyy-mm-dd HH:MM:SS")
    (; hours = hours_arr[:], epoch, lat = lat_arr[:], lon = lon_arr[:])
end

# `RasterStack(url; source=Zarrsource())`'s consolidated-metadata discovery
# only surfaces a subset of this store's ~280 arrays, so named variables are
# opened directly by their store subpath instead.
MicroclimateMapper._contiguous_series_open(url::AbstractString, long_name::AbstractString) =
    Zarr.zopen(url * "/" * long_name)

end
