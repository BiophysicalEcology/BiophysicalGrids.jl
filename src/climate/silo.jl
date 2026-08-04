# SILO (Scientific Information for Land Owners) bindings — daily ~5 km
# gridded weather over Australia, 1889-present. One annual NetCDF per
# variable per year, downloaded and cached via `RasterDataSources.SILO`.
#
# SILO has no native wind field — falls back to CRUCL2's monthly wind
# climatology (`:wind_speed`), same as NicheMapR's own SILO/AWAP implementations.
#
# `:rh_tmax` (RH at the time of Tmax, the day's lowest) maps to
# `reference_humidity_min`; `:rh_tmin` (RH at Tmin, the day's highest) maps
# to `reference_humidity_max` — same convention as GRIDMET's `:rmin`/`:rmax`.

weather_calendar(::Type{<:SILO}) = Daily()
loader(::Type{<:SILO}) = YearlyTimeSeries()
fallback_source(::Type{<:SILO}) = CRUCL2
fallback_layers(::Type{<:SILO}) = (:wind_speed,)
# SILO's DataDrill API requires an email address.
_extra_getpoint_kwargs(::Type{<:SILO}) = (; username = get(ENV, "SILO_EMAIL", ""))

function variables(::Type{<:SILO})
    (
        Variable(Temperature(Maximum()), :max_temp, u"°C"),
        Variable(Temperature(Minimum()), :min_temp, u"°C"),
        # `:daily_rain` is mm depth; numerically equal to kg/m².
        Variable(Rainfall(), :daily_rain, u"kg/m^2"),
        Variable(Reference(RelativeHumidity(Minimum())), :rh_tmax, 1, percent_to_fraction),
        Variable(Reference(RelativeHumidity(Maximum())), :rh_tmin, 1, percent_to_fraction),
        # Daily total solar exposure (MJ/(m²·day))
        Variable(GlobalRadiation(), :radiation, u"MJ/m^2/d"),
    )
end
