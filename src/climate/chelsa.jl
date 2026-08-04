# CHELSA bindings: `CHELSA{Climate}` (1981-2010 monthly climatology) and
# `CHELSA{Future{Climate, ...}}`

# ---------------------------------------------------------------------------
# CHELSA{Climate} — historical baseline climatology
# ---------------------------------------------------------------------------

loader(::Type{<:CHELSA{Climate}}) = MonthlyClimatology()

function variables(::Type{<:CHELSA{Climate}})
    (
        Variable(Temperature(Maximum()), :tasmax, u"°C"),
        Variable(Temperature(Minimum()), :tasmin, u"°C"),
        Variable(WindSpeed(), :sfcWind, u"m/s"),
        Variable(GlobalRadiation(), :rsds, u"W/m^2"),
        Variable(Rainfall(), :pr, u"kg/m^2"),
        # CHELSA gives a single mean relative humidity — use it for both
        # reference_humidity min and max so the vapour-pressure-deficit-based
        # derivation is skipped on both ends.
        Variable(Reference(RelativeHumidity(Minimum())), :hurs, 1, percent_to_fraction),
        Variable(Reference(RelativeHumidity(Maximum())), :hurs, 1, percent_to_fraction),
        # `cloud_cover` is native, so the radiation-based cloud derivation
        # is skipped; cloud_cover_min/max derivations still fire from this value.
        Variable(CloudCover(), :clt, 1, percent_to_fraction),
    )
end

# ---------------------------------------------------------------------------
# CHELSA{Future{Climate, ...}} — projected period climatology
# ---------------------------------------------------------------------------

loader(::Type{<:CHELSA{<:Future{Climate}}}) = MonthlyClimatologyPeriod()

fallback_source(::Type{<:CHELSA{<:Future{Climate}}}) = CHELSA{Climate}

function variables(::Type{<:CHELSA{<:Future{Climate}}})
    (
        Variable(Temperature(Maximum()), :tmax, u"°C"),
        Variable(Temperature(Minimum()), :tmin, u"°C"),
        Variable(Rainfall(), :prec, u"kg/m^2"),
    )
end
