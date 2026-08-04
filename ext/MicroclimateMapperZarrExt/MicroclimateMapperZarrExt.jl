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

end
