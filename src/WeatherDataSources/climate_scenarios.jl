"""
    apply_climate_scenario(Scenario, weather, lon, lat; ystart, yfinish=ystart, kwargs...)

Apply climate change scenario deltas to a weather `NamedTuple` returned by `get_weather`.
Dispatches on the scenario type. `Historical` is a source-agnostic no-op; scenario types
that include a data source (e.g. `TerraClimate{Plus2C}`) control both which deltas to
download and how to apply them.

# Methods
- `apply_climate_scenario(Historical, ...)` — returns `weather` unchanged (source-agnostic)
- `apply_climate_scenario(TerraClimate{Plus2C}, ...)` — applies +2 °C TerraClimate scenario
- `apply_climate_scenario(TerraClimate{Plus4C}, ...)` — applies +4 °C TerraClimate scenario

`weather` may come from any source (TerraClimate, ERA5, etc.) — the source type controls
where the scenario deltas are drawn from, not where the baseline weather came from.

# Variables modified (TerraClimate methods)
- **Temperature** (min/max): additive delta in K
- **Cloud cover** (min/max): scaled by `1 + (1 − srad_ratio)` (inverse solar relationship)
- **Rainfall**: multiplied by precipitation ratio
- **Deep soil temperature**: shifted by annual mean of temperature delta
- **Relative humidity**: actual vapour pressure adjusted for VPD change, then recalculated
  against saturation VP at the new temperature
- **Wind speed**: unchanged

# Example
```julia
weather    = get_weather(TerraClimate, Point([lon, lat]); tstart = Date(2000), elevation)
weather_2c = apply_climate_scenario(TerraClimate{Plus2C}, weather, Point([lon, lat]); tstart = Date(2000))
result_2c  = simulate_microclimate(solar_terrain, micro_terrain, soil_thermal_model, weather_2c)
```
"""
function apply_climate_scenario end

# Historical: no-op regardless of scenario source
function apply_climate_scenario(
    ::Type{Historical},
    weather,
    ::Any;
    vapour_pressure_method = GoffGratch(),
)
    return weather
end

function apply_climate_scenario(
    ::Type{TerraClimate{S}},
    weather,
    point;
    tstart::Date,
    tend::Date = tstart,
    vapour_pressure_method = GoffGratch(),
) where {S <: RasterDataSources.WarmingScenario}
    base = _fetch_tc_raw(point, tstart, tend, Historical)
    scen = _fetch_tc_raw(point, tstart, tend, S)

    # Deltas: temperature additive (K), others multiplicative or additive in native units
    tmax_delta_K = u"K".(scen.tmax .* u"°C") .- u"K".(base.tmax .* u"°C")
    tmin_delta_K = u"K".(scen.tmin .* u"°C") .- u"K".(base.tmin .* u"°C")
    ppt_ratio    = scen.ppt  ./ max.(base.ppt,  1e-6)
    srad_ratio   = scen.srad ./ max.(base.srad, 1e-6)
    vpd_delta    = (scen.vpd .- base.vpd) .* u"kPa"   # positive = drier

    # Reduce to 12-month climatological mean and tile to match weather length
    n_months   = length(weather.environment_minmax.reference_temperature_min)
    tmax_Δ     = _tile_monthly_delta(tmax_delta_K, n_months)
    tmin_Δ     = _tile_monthly_delta(tmin_delta_K, n_months)
    ppt_ratio  = _tile_monthly_delta(ppt_ratio,    n_months)
    srad_ratio = _tile_monthly_delta(srad_ratio,   n_months)
    vpd_Δ      = _tile_monthly_delta(vpd_delta,    n_months)

    em = weather.environment_minmax

    tmin_new = em.reference_temperature_min .+ tmin_Δ
    tmax_new = em.reference_temperature_max .+ tmax_Δ

    rh_max_new = _adjust_rh(em.reference_humidity_max, em.reference_temperature_min,
                             tmin_new, vpd_Δ, vapour_pressure_method)
    rh_min_new = _adjust_rh(em.reference_humidity_min, em.reference_temperature_max,
                             tmax_new, vpd_Δ, vapour_pressure_method)

    cloud_factor  = clamp.(1.0 .+ (1.0 .- srad_ratio), 0.0, 2.0)
    cloud_min_new = clamp.(em.cloud_min .* cloud_factor, 0.0, 1.0)
    cloud_max_new = clamp.(em.cloud_max .* cloud_factor, 0.0, 1.0)

    new_em = MonthlyMinMaxEnvironment(;
        reference_temperature_min = tmin_new,
        reference_temperature_max = tmax_new,
        reference_wind_min        = em.reference_wind_min,
        reference_wind_max        = em.reference_wind_max,
        reference_humidity_min    = clamp.(rh_min_new, 0.0, 1.0),
        reference_humidity_max    = clamp.(rh_max_new, 0.0, 1.0),
        cloud_min                 = cloud_min_new,
        cloud_max                 = cloud_max_new,
        minima_times              = em.minima_times,
        maxima_times              = em.maxima_times,
    )

    ed = weather.environment_daily
    rainfall_new = ed.rainfall .* ppt_ratio
    tmean_Δ      = (tmax_Δ .+ tmin_Δ) ./ 2
    nyears_w     = max(1, length(ed.deep_soil_temperature) ÷ 12)
    deep_T_new   = ed.deep_soil_temperature .+ _annual_means_monthly(tmean_Δ, nyears_w)

    new_ed = DailyTimeseries(;
        shade              = ed.shade,
        soil_wetness       = ed.soil_wetness,
        surface_emissivity = ed.surface_emissivity,
        cloud_emissivity   = ed.cloud_emissivity,
        rainfall           = rainfall_new,
        deep_soil_temperature = deep_T_new,
        leaf_area_index    = ed.leaf_area_index,
    )

    return merge(weather, (; environment_minmax = new_em, environment_daily = new_ed))
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _fetch_tc_raw(point, tstart::Date, tend::Date, scenario)
    layers = (:tmax, :tmin, :ppt, :vpd, :srad)
    tmax = Float64[]; tmin = Float64[]
    ppt  = Float64[]; vpd  = Float64[]; srad = Float64[]
    for yr in year(tstart):year(tend)
        paths = getraster(TerraClimate{scenario}, layers; date = Date(yr))
        for (sym, vec) in ((:tmax, tmax), (:tmin, tmin),
                           (:ppt, ppt), (:vpd, vpd), (:srad, srad))
            append!(vec, _extract_monthly(Raster(paths[sym]; lazy = true), point))
        end
    end
    return (; tmax, tmin, ppt, vpd, srad)
end

function _tile_monthly_delta(delta::AbstractVector, n::Int)
    mean_monthly = [mean(delta[m:12:end]) for m in 1:12]
    return repeat(mean_monthly, cld(n, 12))[1:n]
end

function _adjust_rh(
    rh_orig::AbstractVector,
    T_orig::AbstractVector,
    T_new::AbstractVector,
    vpd_delta_kPa::AbstractVector,
    method,
)
    return map(zip(rh_orig, T_orig, T_new, vpd_delta_kPa)) do (rh, T_o, T_n, dvpd)
        e_sat_orig_kPa   = ustrip(u"hPa", vapour_pressure(method, T_o)) / 10.0
        e_actual_kPa     = rh * e_sat_orig_kPa
        e_actual_new_kPa = max(e_actual_kPa - ustrip(u"kPa", dvpd), 0.0)
        e_sat_new_kPa    = ustrip(u"hPa", vapour_pressure(method, T_n)) / 10.0
        e_sat_new_kPa <= 0.0 && return 0.0
        clamp(e_actual_new_kPa / e_sat_new_kPa, 0.0, 1.0)
    end
end
