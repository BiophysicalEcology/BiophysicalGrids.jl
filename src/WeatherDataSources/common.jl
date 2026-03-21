const MONTHLY_BASE_DAYS = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]

"""
    monthly_days(nyears) → Vector{Int}

Day-of-year values for the midpoint of each month, repeated for `nyears` years.
Used to set the `days` field of `MicroProblem` for monthly forcing data.
"""
monthly_days(nyears::Int) = repeat(MONTHLY_BASE_DAYS, nyears)

"""
    mean_annual_temperature(tmax_vec, tmin_vec) → Quantity

Compute annual mean temperature from vectors of monthly max and min temperatures.
Used as an estimate of deep soil temperature.
"""
function mean_annual_temperature(tmax_vec::AbstractVector, tmin_vec::AbstractVector)
    return mean((tmax_vec .+ tmin_vec) ./ 2)
end

"""
    extract_point(raster, lon, lat)

Extract a scalar (or time-series) value from `raster` at the nearest grid point
to (`lon`, `lat`). Returns the bare array when the raster has a time dimension.
"""
function extract_point(raster, lon::Real, lat::Real)
    return raster[X(Near(lon)), Y(Near(lat))]
end

# ---------------------------------------------------------------------------
# Generic simulation entry point
# ---------------------------------------------------------------------------

"""
    simulate_microclimate(solar_terrain, micro_terrain, soil_thermal_model, weather; kwargs...)

Run the microclimate model for a single location.

`weather` is a `NamedTuple` returned by `get_weather`, containing at minimum:
- `environment_minmax::MonthlyMinMaxEnvironment`
- `environment_daily::DailyTimeseries`
- `environment_hourly::HourlyTimeseries` (or `nothing`)
- `latitude` — latitude with Unitful degrees
- `days` — day-of-year vector matching the environment data length

The user is responsible for constructing the terrain and soil model structs with
appropriate parameters for the simulation site.

# Keyword arguments
- `soil_moisture_model`: if `nothing` (default), a sensible default is built from
  `soil_thermal_model.bulk_density` and `soil_thermal_model.mineral_density`.
- `depths`: soil node depths (default: Microclimate.jl's 19-node `DEFAULT_DEPTHS`).
- `heights`: air profile node heights (default: `[0.01, 2.0]u"m"`).
- `runmoist`: enable soil moisture simulation (default: `false`).
- `solar_model`: `SolarProblem` instance (default: `SolarProblem()`).
- `iterate_day`: number of iterations per day (default: 3).
- `spinup`: spin up the first day (default: `false`).
- `initial_soil_temperature`: initial soil temperatures. Defaults to `nothing`,
  which lets Microclimate.jl use the mean air temperature of each month as the
  starting T₀ — matching NicheMapR's monthly (`microdaily=0`) behaviour.
- `initial_soil_moisture`: initial volumetric soil moisture fractions.

# Returns
`MicroResult` from `Microclimate.solve`.
"""
function simulate_microclimate(
    solar_terrain::SolarTerrain,
    micro_terrain::MicroTerrain,
    soil_thermal_model,
    weather::NamedTuple;
    soil_moisture_model = nothing,
    depths = Microclimate.DEFAULT_DEPTHS,
    heights = [0.01, 2.0]u"m",
    runmoist::Bool = false,
    solar_model::SolarProblem = SolarProblem(),
    iterate_day::Int = 3,
    spinup::Bool = false,
    initial_soil_temperature = nothing,
    initial_soil_moisture = fill(0.42 * 0.25, length(depths)),
    kwargs...,
)
    (; environment_minmax, environment_daily, environment_hourly, latitude, days) = weather
    # Build (ndepths × ndays) precomputed soil moisture matrix from monthly weather data.
    # Used by Microclimate.jl when runmoist=false to vary soil moisture per day.
    precomputed_soil_moisture = let sm = get(weather, :soil_moisture_monthly, nothing)
        if isnothing(sm)
            nothing
        else
            # Repeat each monthly value for all depth nodes (uniform with depth).
            # NicheMapR only sets top 2 nodes; using uniform is more physically consistent.
            repeat(sm', length(depths), 1)  # (ndepths × nmonths)
        end
    end

    # Build default soil moisture model from soil thermal parameters if not provided
    if isnothing(soil_moisture_model)
        soil_moisture_model = example_soil_moisture_model(
            depths;
            bulk_density = ustrip(u"Mg/m^3", soil_thermal_model.bulk_density),
            mineral_density = ustrip(u"Mg/m^3", soil_thermal_model.mineral_density),
            root_density = fill(0.0, length(depths))u"m/m^3",
        )
    end

    # initial_soil_temperature = nothing → Microclimate.jl resets T0 to the mean
    # air temperature of each month (NicheMapR microdaily=0 behaviour).

    problem = MicroProblem(;
        latitude,
        days,
        hours = collect(0.0:1.0:23.0),
        depths,
        heights,
        solar_model,
        solar_terrain,
        micro_terrain,
        soil_moisture_model,
        soil_thermal_model,
        environment_minmax,
        environment_daily,
        environment_hourly,
        iterate_day,
        runmoist,
        spinup,
        initial_soil_temperature,
        initial_soil_moisture,
        precomputed_soil_moisture,
        kwargs...,
    )

    return Microclimate.solve(problem)
end
