# ERA5 bindings — hourly global reanalysis from ARCO-ERA5

weather_calendar(::Type{<:ERA5}) = Daily()
native_timestep(::Type{<:ERA5}) = Hourly()
loader(::Type{<:ERA5}) = ContiguousTimeSeries()

function variables(::Type{<:ERA5})
    (
        Variable(Reference(Temperature()), :t2m, u"K"),
        Variable(EastwardWindSpeed(), :u10, u"m/s"),
        Variable(NorthwardWindSpeed(), :v10, u"m/s"),
        Variable(DewpointTemperature(), :d2m, u"K"),
        Variable(Pressure(), :sp, u"Pa"),
        # Total cloud cover already 0-1 fraction; no unit, no transform.
        Variable(CloudCover(), :tcc, 1),
        Variable(GlobalRadiation(), :ssrd, u"J/m^2/hr"),
        Variable(LongwaveRadiation(), :strd, u"J/m^2/hr"),
        # `:tp` is total precipitation in metres of water per hour; multiply
        # by 1000 to get kg/m² (assuming water density of 1000 kg/m³).
        Variable(Rainfall(), :tp, u"kg/m^2", raw -> raw * 1000.0),
    )
end

# Not part of `variables()` (would pull the full requested Ti range at
# every timestep) -- only used by `_load_init_snapshot` for init seeding.
init_variables(::Type{<:ERA5}) = (
    variables(ERA5)...,
    Variable(SoilTemperature(Mean()), :stl1, u"K"),
    Variable(SoilMoisture(), :swvl1, 1),
)
