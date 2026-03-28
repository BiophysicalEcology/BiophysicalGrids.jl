"""
    get_weather(ERA5, point; tstart, tend, kwargs...) -> NamedTuple

Download ERA5 hourly reanalysis data for a point location via the ARCO-ERA5 cloud
Zarr store (requires an internet connection; data are cached locally by `Zarr.jl`).

All `HourlyTimeseries` fields are populated directly from ERA5 hourly values.
Returns a `DailyMinMaxEnvironment` (per-day min/max derived from the hourly data)
which triggers `daily=true` mode in `simulate_microclimate` — consecutive calendar
days inherit soil state and iterate once rather than converging independently.

# Arguments
- `point` : GeoInterface-compatible point geometry (e.g. `Point([lon, lat])`);
  longitude in −180..360, latitude in −90..90
- `tstart`, `tend` : `DateTime` bounds (inclusive; UTC; must span whole days)

# Keyword arguments
- `elevation` : site elevation as a Unitful quantity (e.g. `270.0u"m"`).  When
  provided together with `grid_elevation`, temperatures are lapse-rate-corrected
  from the ERA5 grid-cell elevation to the site elevation.
- `grid_elevation` : ERA5 model orography at the nearest grid cell (m).  Defaults
  to `elevation` (no correction) when not supplied.
- `lapse_rate_type` : `LapseRate` subtype, default `EnvironmentalLapseRate()`.
- `vapour_pressure_method` : FluidProperties vapour pressure equation used to
  convert ERA5 dewpoint → RH (default `GoffGratch()`).
- `use_era5_soil_moisture` : if `true`, extract `swvl1` (0–7 cm volumetric soil
  water, m³/m³) and use as soil moisture forcing.  Default `false`.

# Returns
`NamedTuple` with fields:
- `environment_minmax::DailyMinMaxEnvironment` — per-day min/max
- `environment_daily::DailyTimeseries`
- `environment_hourly::HourlyTimeseries` — all fields populated
- `latitude` — site latitude with `u"°"` units
- `days` — day-of-year integers for each simulated day
- `soil_moisture_monthly` — daily soil moisture vector (m³/m³), uniform with depth

# Example
```julia
using BiophysicalGrids, Dates
weather = get_weather(ERA5, Point([-89.46, 43.14]);
    tstart    = DateTime(2000, 1, 1),
    tend      = DateTime(2000, 12, 31, 23),
    elevation = 270.0u"m",
)
```
"""
function get_weather(
    ::Type{ERA5},
    point;
    tstart::DateTime,
    tend::DateTime,
    elevation = nothing,
    grid_elevation = nothing,
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
    vapour_pressure_method = GoffGratch(),
    use_era5_soil_moisture::Bool = false,
)
    lon, lat = _lonlat(point)
    # -----------------------------------------------------------------------
    # Open ERA5 ARCO Zarr store
    # -----------------------------------------------------------------------
    source = getraster(ERA5)
    store  = Zarr.CachingHTTPStore(source)
    ds     = zopen(Zarr.ConsolidatedStore(store, ""))

    # -----------------------------------------------------------------------
    # Spatial indices
    # -----------------------------------------------------------------------
    lats    = ds["latitude"][:]
    lons    = ds["longitude"][:]
    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- mod(lon, 360.0)))

    # -----------------------------------------------------------------------
    # Time index — compute directly from ERA5 epoch (hours since 1900-01-01)
    # to avoid reading the full time axis (~750k values for the complete archive).
    # ERA5 is hourly and contiguous, so the index is exact arithmetic.
    # -----------------------------------------------------------------------
    hours_since_1900(dt) = Dates.value(dt - DateTime(1900, 1, 1)) ÷ 3_600_000
    t_first     = Int(first(ds["time"][1:1]))  # single-chunk read: first archived hour
    t_start_idx = hours_since_1900(tstart) - t_first + 1
    t_end_idx   = hours_since_1900(tend)   - t_first + 1
    t_start_idx < 1 &&
        error("tstart=$tstart is before the ERA5 archive start")
    t_end_idx > size(ds["time"], 1) &&
        error("tend=$tend is after the ERA5 archive end")

    t_idx  = t_start_idx:t_end_idx
    nhours = length(t_idx)
    ndays  = nhours ÷ 24
    ndays == 0 && error("ERA5 time range covers fewer than 24 hours")
    nhours != ndays * 24 &&
        @warn "ERA5 range spans $nhours hours; truncating to $ndays complete days"
    use_hours = 1:ndays*24   # only whole days

    # -----------------------------------------------------------------------
    # Extract ERA5 variables — stored as (lon × lat × time)
    # -----------------------------------------------------------------------
    extr(sym) = Float64.(vec(ds[layername(ERA5, sym)][lon_idx:lon_idx, lat_idx:lat_idx, t_idx]))[use_hours]

    t2m  = extr(:t2m)   # K
    d2m  = extr(:d2m)   # K  (dewpoint at 2 m)
    u10  = extr(:u10)   # m/s
    v10  = extr(:v10)   # m/s
    sp   = extr(:sp)    # Pa
    ssrd = extr(:ssrd)  # J/m²/hr (accumulated per hour)
    strd = extr(:strd)  # J/m²/hr (accumulated per hour)
    tcc  = extr(:tcc)   # 0–1 fraction
    tp   = extr(:tp)    # m/hr

    # -----------------------------------------------------------------------
    # Lapse-rate correction on 2m temperature
    # -----------------------------------------------------------------------
    site_elev = isnothing(elevation) ? 0.0u"m" : elevation
    ref_elev  = isnothing(grid_elevation) ? site_elev : grid_elevation
    Δz        = site_elev - ref_elev

    t2m_K = t2m .* u"K"
    if !iszero(Δz)
        t2m_K = lapse_adjust_temperature(t2m_K, Δz, lapse_rate_type)
    end

    # -----------------------------------------------------------------------
    # Dewpoint → relative humidity
    # -----------------------------------------------------------------------
    _vp(T) = ustrip(u"hPa", vapour_pressure(vapour_pressure_method, T))
    d2m_K  = d2m .* u"K"
    rh_hr  = clamp.(_vp.(d2m_K) ./ _vp.(t2m_K), 0.0, 1.0)

    # -----------------------------------------------------------------------
    # Wind speed: 10 m components → 2 m scalar (power-law height correction)
    # -----------------------------------------------------------------------
    wind_hr = @. sqrt(u10^2 + v10^2) * (2.0 / 10.0)^0.15 * u"m/s"

    # -----------------------------------------------------------------------
    # Radiation: J/m²/hr accumulated → W/m²  (clamp to ≥ 0)
    # -----------------------------------------------------------------------
    solar_hr = max.(ssrd, 0.0) ./ 3600.0 .* u"W/m^2"
    lw_hr    = max.(strd, 0.0) ./ 3600.0 .* u"W/m^2"

    # -----------------------------------------------------------------------
    # Precipitation: m/hr → kg/m²/hr  (clamp to ≥ 0)
    # -----------------------------------------------------------------------
    rain_hr = max.(tp, 0.0) .* 1000.0 .* u"kg/m^2"

    # -----------------------------------------------------------------------
    # Daily aggregations via reshape (24 hours × ndays)
    # -----------------------------------------------------------------------
    function daily_ext(v, f)
        mat = reshape(v, 24, ndays)
        return [f(mat[:, d]) for d in 1:ndays]
    end
    function daily_sum(v)
        mat = reshape(v, 24, ndays)
        return vec(sum(mat; dims = 1))
    end

    # Day-of-year vector
    dates = [Date(tstart) + Day(d - 1) for d in 1:ndays]
    doys  = dayofyear.(dates)

    # Deep soil temperature: mean of all daily mean temperatures
    daily_tmean = [mean(t2m_K[((d-1)*24+1):d*24]) for d in 1:ndays]
    deep_soil_T = fill(mean(daily_tmean), ndays)

    # Daily rainfall totals (sum of 24 hourly values)
    daily_rain = daily_sum(ustrip.(u"kg/m^2", rain_hr)) .* u"kg/m^2"

    # -----------------------------------------------------------------------
    # Optional ERA5 soil moisture (swvl1: 0–7 cm, m³/m³)
    # -----------------------------------------------------------------------
    soil_moisture_monthly = if use_era5_soil_moisture
        swvl1 = extr(:swvl1)
        [mean(swvl1[((d-1)*24+1):d*24]) for d in 1:ndays]
    else
        fill(0.42 * 0.25, ndays)
    end

    # -----------------------------------------------------------------------
    # Build Microclimate.jl environment structs
    # -----------------------------------------------------------------------
    wind_ms = ustrip.(u"m/s", wind_hr)

    environment_minmax = DailyMinMaxEnvironment(;
        reference_temperature_min = daily_ext(t2m_K, minimum),
        reference_temperature_max = daily_ext(t2m_K, maximum),
        reference_wind_min        = daily_ext(wind_ms, minimum) .* u"m/s",
        reference_wind_max        = daily_ext(wind_ms, maximum) .* u"m/s",
        reference_humidity_min    = daily_ext(rh_hr, minimum),
        reference_humidity_max    = daily_ext(rh_hr, maximum),
        cloud_min                 = daily_ext(tcc, minimum),
        cloud_max                 = daily_ext(tcc, maximum),
        minima_times = (temp = 0, wind = 0, humidity = 1, cloud = 1),
        maxima_times = (temp = 1, wind = 1, humidity = 0, cloud = 0),
    )

    environment_daily = DailyTimeseries(;
        shade              = fill(0.0, ndays),
        soil_wetness       = fill(0.0, ndays),
        surface_emissivity = fill(0.95, ndays),
        cloud_emissivity   = fill(0.95, ndays),
        rainfall           = daily_rain,
        deep_soil_temperature = deep_soil_T,
        leaf_area_index    = fill(0.1, ndays),
    )

    environment_hourly = HourlyTimeseries(;
        pressure              = sp .* u"Pa",
        reference_temperature = t2m_K,
        reference_humidity    = rh_hr,
        reference_wind_speed  = wind_hr,
        global_radiation      = solar_hr,
        longwave_radiation    = lw_hr,
        cloud_cover           = tcc,
        rainfall              = rain_hr,
        zenith_angle          = nothing,
    )

    return (;
        environment_minmax,
        environment_daily,
        environment_hourly,
        latitude = lat * u"°",
        days     = doys,
        soil_moisture_monthly,
    )
end
