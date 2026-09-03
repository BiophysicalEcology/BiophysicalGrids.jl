# BARRA bindings — Bureau of Meteorology's high-resolution Australian
# regional reanalysis (BARRA-R2 ~11 km / BARRA-C2 ~4 km), Australia's
# analog to ERA5. One NetCDF per variable per real calendar month.

weather_calendar(::Type{<:BARRA}) = Daily()
native_timestep(::Type{<:BARRA}) = Hourly()
loader(::Type{<:BARRA}) = MonthlyTimeSeries()

function variables(::Type{<:BARRA})
    (
        Variable(Reference(Temperature()),      :tas,     u"K"),
        Variable(WindSpeed(),                   :sfcWind, u"m/s"),
        Variable(Reference(RelativeHumidity()), :hurs,    1, percent_to_fraction),
        Variable(SeaLevelPressure(),            :psl,     u"Pa"),
        Variable(GlobalRadiation(),             :rsds,    u"W/m^2"),
        Variable(LongwaveRadiation(),           :rlds,    u"W/m^2"),
        # `:pr` (kg/m²/s) × 3600 s → kg/m² per hour.
        Variable(Rainfall(),                    :pr,      u"kg/m^2", raw -> raw * 3600.0),
        Variable(Elevation(),                   :orog,    u"m"),
    )
end

# BARRA's own `:orog` grid as the DEM — same role as CRUCL2's `:elv`. 0.3°
# buffer matches SRTM's, comfortably wider than a BARRA-R2 cell (~11 km).
# `:orog` only exists under the `Static` frequency parameter (RasterDataSources
# `BARRA{P,D,Static}`), not under `T`'s own (possibly Hourly/Daily/Monthly) one.
function _load_dem(T::Type{<:BARRA{P,D}}, area::Extent) where {P,D}
    buffered = Extents.buffer(area, (X = 0.3, Y = 0.3))
    path = getraster(BARRA{P,D,RasterDataSources.Static}, :orog)
    orog = read(crop(Raster(path; name = :orog, lazy = true); to = buffered, touches = true))
    return Rasters.replace_missing(orog, 0)
end
