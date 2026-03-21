"""
    get_weather(TerraClimate, lon, lat; ystart, yfinish, kwargs...)
    get_weather(TerraClimate{Plus2C}, lon, lat; ...)
    get_weather(TerraClimate{Plus4C}, lon, lat; ...)

Download TerraClimate monthly climate data for the given location and time period and
return a `NamedTuple` of Microclimate.jl environment objects ready to pass to
`simulate_microclimate`.

TerraClimate provides global monthly data at ~4 km resolution (1958–2024 for historical).
Variables extracted: `tmax`, `tmin`, `ws` (wind), `vap` (vapour pressure), `vpd` (vapour
pressure deficit), `srad` (shortwave radiation), `ppt` (precipitation).

# Arguments
- `lon`, `lat`: decimal degrees (longitude, latitude)
- `ystart`, `yfinish`: simulation years (inclusive). Default `yfinish = ystart`.

# Keyword arguments
- `elevation`: site elevation as a Unitful quantity (e.g. `270.0u"m"`). If provided,
  temperatures are lapse-rate-corrected from `grid_elevation` to `elevation`, and
  atmospheric pressure is computed at site elevation. If `nothing`, no lapse correction.
- `grid_elevation`: TerraClimate grid-cell elevation at `(lon, lat)`, used as the
  reference for lapse rate correction. Defaults to `elevation` (no correction) when
  not provided — NicheMapR's equivalent of `elev = NA`. For accurate adjustment supply
  the WorldClim 2.5-arcmin elevation at the grid cell (the same reference NicheMapR
  uses internally), or an SRTM-derived value from `Geomorphometry`.
- `lapse_rate_type::LapseRate`: temperature correction type
  (default: `EnvironmentalLapseRate()`).
- `vapour_pressure_method`: FluidProperties vapour pressure equation for
  `vap`/`vpd` → relative humidity conversion (default: `GoffGratch()`).
  Alternatives: `Teten()`, `Huang()`.

# Returns
`NamedTuple` with fields:
- `environment_minmax::MonthlyMinMaxEnvironment`
- `environment_daily::DailyTimeseries`
- `environment_hourly::HourlyTimeseries`
- `latitude` — site latitude with `u"°"` units
- `days` — day-of-year midpoints repeated for `nyears`

# Example
```julia
weather = get_weather(TerraClimate, -89.4557, 43.1379;
    ystart = 2000, yfinish = 2001,
    elevation = 270.0u"m",
    vapour_pressure_method = GoffGratch(),
    lapse_rate_type = EnvironmentalLapseRate(),
)
```
"""
function get_weather(
    ::Type{<:TerraClimate},
    lon::Real,
    lat::Real;
    ystart::Int,
    yfinish::Int = ystart,
    scenario::Type{<:RasterDataSources.WarmingScenario} = Historical,
    elevation = nothing,
    grid_elevation = nothing,
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
    vapour_pressure_method = GoffGratch(),
)
    nyears = yfinish - ystart + 1
    years  = ystart:yfinish

    layers = (:tmax, :tmin, :ws, :vap, :vpd, :srad, :ppt, :soil)

    # Accumulate monthly values across years
    tmax_all = Float64[]
    tmin_all = Float64[]
    ws_all   = Float64[]
    vap_all  = Float64[]
    vpd_all  = Float64[]
    srad_all = Float64[]
    ppt_all  = Float64[]
    soil_all = Float64[]

    for yr in years
        paths = getraster(TerraClimate{scenario}, layers; date = Date(yr))
        for (sym, vec) in (
            (:tmax, tmax_all), (:tmin, tmin_all), (:ws, ws_all),
            (:vap, vap_all), (:vpd, vpd_all), (:srad, srad_all), (:ppt, ppt_all),
            (:soil, soil_all),
        )
            r    = Raster(paths[sym]; lazy = true)
            vals = _extract_monthly(r, lon, lat)
            append!(vec, vals)
        end
    end

    nmonths = length(tmax_all)  # should be 12 * nyears

    # Elevations and atmospheric pressure
    site_elev = isnothing(elevation) ? 0.0u"m" : elevation
    # Default grid_elevation to site_elev (no lapse correction) when not provided.
    # NicheMapR uses WorldClim DEM elevation as the TC grid reference; without it,
    # assuming sea level (0m) over-corrects most sites. Provide grid_elevation
    # explicitly (e.g. from WorldClim or SRTM) for accurate lapse adjustment.
    ref_elev  = isnothing(grid_elevation) ? site_elev : grid_elevation
    Δz        = site_elev - ref_elev
    P_atm     = atmospheric_pressure(site_elev)

    # Temperature with optional lapse rate correction.
    # Convert to K immediately so arithmetic (mean, addition) works regardless of
    # whether lapse correction is applied — affine °C cannot be summed in Unitful.
    tmax_C = u"K".(tmax_all .* u"°C")
    tmin_C = u"K".(tmin_all .* u"°C")
    if !iszero(Δz)
        tmax_C = lapse_adjust_temperature(tmax_C, Δz, lapse_rate_type)
        tmin_C = lapse_adjust_temperature(tmin_C, Δz, lapse_rate_type)
    end

    # Relative humidity — match NicheMapR micro_terra.R approach:
    # actual VP = e_sat(Tmean) − VPD  (VPD defined at mean temperature)
    # RH_max = actual_VP / e_sat(Tmin)   (highest humidity at coldest time)
    # RH_min = actual_VP / e_sat(Tmax)   (lowest humidity at warmest time)
    tmean_C = (tmax_C .+ tmin_C) ./ 2
    rh_max = _rh_from_vpd_at_tmean(vpd_all, tmean_C, tmin_C, vapour_pressure_method)
    rh_min = _rh_from_vpd_at_tmean(vpd_all, tmean_C, tmax_C, vapour_pressure_method)

    # Wind — match NicheMapR micro_terra.R:
    # TerraClimate wind is at 10 m; correct to 2 m using power-law wind shear
    # (exponent 0.15 for open/grassland terrain per NicheMapR default).
    # Diurnal cycle: max = corrected mean, min = 10% of max.
    _wind_h_factor = (2.0 / 10.0)^0.15   # ≈ 0.794
    ws_max = ws_all .* _wind_h_factor .* u"m/s"
    ws_min = ws_max .* 0.1

    # Cloud cover from shortwave radiation vs. clear-sky estimate
    # Use a flat-terrain SolarTerrain at the site for the clear-sky baseline
    _solar_terrain_flat = SolarTerrain(;
        elevation            = site_elev,
        slope                = 0.0u"°",
        aspect               = 0.0u"°",
        horizon_angles       = fill(0.0u"°", 24),
        albedo               = 0.15,
        atmospheric_pressure = P_atm,
        latitude             = lat * u"°",
        longitude            = lon * u"°",
    )
    doys       = repeat(MONTHLY_BASE_DAYS, nyears)
    cloud_vals = cloud_from_srad(srad_all .* u"W/m^2", _solar_terrain_flat, doys)

    # Deep soil temperature: annual mean for each year, broadcast to monthly
    tmid_all   = (tmax_C .+ tmin_C) ./ 2
    deep_soil_T = _annual_means_monthly(tmid_all, nyears)

    # Precipitation: TerraClimate ppt in mm/month; mm ≡ kg/m²
    rainfall = ppt_all .* u"kg/m^2"

    # Soil moisture: TerraClimate `soil` = total column (0–2 m) water in mm.
    # Convert to volumetric water content (m³/m³) using NicheMapR formula:
    #   θ = soil_mm / 1000 × (1 − bulk_density / mineral_density)
    # Defaults match NicheMapR's BD=1.3 Mg/m³, MD=2.56 Mg/m³ → porosity ≈ 0.492.
    # Returns a (nmonths,) vector; broadcast to (ndepths × nmonths) in simulate_microclimate.
    soil_moisture_monthly = soil_all ./ 1000.0 .* (1.0 - 1.3 / 2.56)  # mm → m³/m³, porosity ≈ 0.492

    # Pressure timeseries: constant per hour × day (nmonths × 24 hourly values)
    pressure_ts = fill(P_atm, nmonths * 24)

    # -----------------------------------------------------------------------
    # Build Microclimate.jl environment structs
    # -----------------------------------------------------------------------
    environment_minmax = MonthlyMinMaxEnvironment(;
        reference_temperature_min = tmin_C,
        reference_temperature_max = tmax_C,
        reference_wind_min        = ws_min,
        reference_wind_max        = ws_max,
        reference_humidity_min    = clamp.(rh_min, 0.0, 1.0),
        reference_humidity_max    = clamp.(rh_max, 0.0, 1.0),
        cloud_min                 = Float64.(cloud_vals),
        cloud_max                 = Float64.(cloud_vals),
        minima_times = (temp = 0, wind = 0, humidity = 1, cloud = 1),
        maxima_times = (temp = 1, wind = 1, humidity = 0, cloud = 0),
    )

    environment_daily = DailyTimeseries(;
        shade              = fill(0.0, nmonths),
        soil_wetness       = fill(0.0, nmonths),
        surface_emissivity = fill(0.95, nmonths),
        cloud_emissivity   = fill(0.95, nmonths),
        rainfall,
        deep_soil_temperature = deep_soil_T,
        leaf_area_index    = fill(0.1, nmonths),
    )

    environment_hourly = HourlyTimeseries(;
        pressure              = pressure_ts,
        reference_temperature = nothing,
        reference_humidity    = nothing,
        reference_wind_speed  = nothing,
        global_radiation      = nothing,
        longwave_radiation    = nothing,
        cloud_cover           = nothing,
        rainfall              = nothing,
        zenith_angle          = nothing,
    )

    return (;
        environment_minmax,
        environment_daily,
        environment_hourly,
        latitude = lat * u"°",
        days     = doys,
        soil_moisture_monthly,   # (nmonths,) volumetric water content m³/m³
    )
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""Extract the 12 monthly values from a TerraClimate raster at the nearest grid point."""
function _extract_monthly(raster, lon::Real, lat::Real)
    ts = raster[X(Near(lon)), Y(Near(lat))]
    return Float64.(collect(ts))
end

"""
    _rh_from_vpd_at_tmean(vpd_kPa, Tmean_vec, Tref_vec, method) → Vector{Float64}

Derive relative humidity matching the NicheMapR micro_terra.R approach:
  actual_VP = e_sat(Tmean) − VPD
  RH        = actual_VP / e_sat(Tref)

Pass `Tref = Tmin` for RH_max and `Tref = Tmax` for RH_min.
FluidProperties `vapour_pressure` returns hPa; converted to kPa here.
"""
function _rh_from_vpd_at_tmean(
    vpd_kPa::AbstractVector,
    Tmean_vec::AbstractVector,
    Tref_vec::AbstractVector,
    method,
)
    return map(zip(vpd_kPa, Tmean_vec, Tref_vec)) do (vpd, Tmean, Tref)
        e_s_mean_kPa = ustrip(u"hPa", vapour_pressure(method, Tmean)) / 10.0
        actual_e_kPa = e_s_mean_kPa - vpd
        actual_e_kPa <= 0.0 && return 0.0
        e_s_ref_kPa  = ustrip(u"hPa", vapour_pressure(method, Tref))  / 10.0
        e_s_ref_kPa  <= 0.0 && return 1.0
        clamp(actual_e_kPa / e_s_ref_kPa, 0.0, 1.0)
    end
end

"""
Broadcast annual mean temperatures to monthly resolution.
For each of `nyears` years (12 months each), every month in the year gets that
year's annual mean temperature. Used for `deep_soil_temperature`.
"""
function _annual_means_monthly(tmid_all::AbstractVector, nyears::Int)
    result = similar(tmid_all)
    for y in 1:nyears
        idx = (1:12) .+ 12 * (y - 1)
        yr_mean = mean(tmid_all[idx])
        result[idx] .= yr_mean
    end
    return result
end
