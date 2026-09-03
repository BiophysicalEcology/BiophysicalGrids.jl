# ECMWFERA5/ECMWFERA5Land bindings — ECMWF's own authenticated ARCO Zarr
# stores (needs a CDS API key, same as CDSERA5), queue-free unlike CDSERA5.

weather_calendar(::Type{<:ECMWFERA5}) = Daily()
native_timestep(::Type{<:ECMWFERA5}) = Hourly()
# Not ContiguousTimeSeries: unlike GCP's monolithic ERA5, getraster(ECMWFERA5)
# (no layer) resolves to RasterDataSources.jl's generic `getraster(T) =
# getraster(T, layers(T))` fallback, and layers(ECMWFERA5) is its topic
# groups (:sfc, :wav) -- a NamedTuple of both, not one usable source.
loader(::Type{<:ECMWFERA5}) = MultiGroupContiguousTimeSeries()
native_group(::Type{<:ECMWFERA5}, ::Symbol) = :sfc  # every field lives here, not :wav
# Default points_load_buffer (2°) is tuned for ~1.9° grids; at ERA5's ~0.25°
# spacing that pulls in an 8x8 cell block per point instead of a handful.
points_load_buffer(::Type{<:ECMWFERA5}) = 0.5

function variables(::Type{<:ECMWFERA5})
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

init_variables(::Type{<:ECMWFERA5}) = variables(ECMWFERA5)

weather_calendar(::Type{<:ECMWFERA5Land}) = Daily()
native_timestep(::Type{<:ECMWFERA5Land}) = Hourly()
loader(::Type{<:ECMWFERA5Land}) = MultiGroupContiguousTimeSeries()
# At ERA5-Land's ~9km (~0.1°) spacing, the 2° default pulls a ~40x40 cell
# block per point -- hundreds of MB of chunks for a single point's series.
points_load_buffer(::Type{<:ECMWFERA5Land}) = 0.2
fallback_source(::Type{<:ECMWFERA5Land}) = ECMWFERA5
fallback_layers(::Type{<:ECMWFERA5Land}) = (:cloud_cover,)  # Land has no cloud-cover group

function variables(::Type{<:ECMWFERA5Land})
    (
        Variable(Reference(Temperature()), :t2m, u"K"),
        Variable(DewpointTemperature(), :d2m, u"K"),
        Variable(EastwardWindSpeed(), :u10, u"m/s"),
        Variable(NorthwardWindSpeed(), :v10, u"m/s"),
        Variable(Pressure(), :sp, u"Pa"),
        Variable(Rainfall(), :tp, u"kg/m^2", raw -> raw * 1000.0),
        Variable(GlobalRadiation(), :ssrd, u"J/m^2/hr"),
        Variable(LongwaveRadiation(), :strd, u"J/m^2/hr"),
    )
end

# Which ECMWFERA5Land Zarr store each field lives in -- all groupings
# confirmed live.
function native_group(::Type{<:ECMWFERA5Land}, f::Symbol)
    f in (:t2m, :d2m) && return :sfc_2m_temperature
    f in (:u10, :v10) && return :sfc_wind
    f in (:sp, :tp) && return :sfc_pressure_precipitation
    f in (:ssrd, :strd) && return :sfc_radiation_heat
    f === :stl1 && return :sfc_soil_temperature
    f === :swvl1 && return :sfc_soil_water
    error("no ECMWFERA5Land group declared for field :$f")
end

init_variables(::Type{<:ECMWFERA5Land}) = (
    variables(ECMWFERA5Land)...,
    Variable(SoilTemperature(Mean()), :stl1, u"K"),
    Variable(SoilMoisture(), :swvl1, 1),
)
