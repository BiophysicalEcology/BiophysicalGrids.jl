# Per-pixel gridded solar radiation calculation.
#
# Provides `solar_radiation_grid`, which runs `SolarRadiation.solar_radiation`
# for every pixel in a terrain grid across a set of days and hours.

"""
    solar_radiation_grid(
        terrain_grids, solar_model, days, hours;
        albedo, extract, verbose
    ) -> Array

Compute solar radiation for every pixel in a terrain grid across the given
`days` and `hours`.

# Positional arguments
- `terrain_grids` : NamedTuple from [`compute_terrain_grids`](@ref)
- `solar_model`   : [`SolarProblem`](@ref) instance
- `days`          : vector of day-of-year values (e.g. `[172.0]`)
- `hours`         : vector of hour values (e.g. `6.0:1.0:20.0`)

# Keyword arguments
- `albedo`   : scalar or `(ny, nx)` matrix (default `0.2`)
- `extract`  : function `(result) -> value` applied to each pixel's
  `SolarResult` to extract the output scalar (default: `r -> r.global_terrain[1]`,
  which returns global terrain radiation in W/m² for the first day/hour).
  For multi-day/hour runs supply a custom extractor that returns a `Vector`.
- `verbose`  : print per-hour progress (default `true`)

# Returns
`Array{Union{Float64,Missing},N}` of shape `(ny, nx)` when `extract` returns a
scalar per pixel, or `(ny, nx, n)` when it returns a length-`n` vector.
Pixels with missing terrain data are `missing`.

# Example — daily total radiation on the summer solstice
```julia
hours_vec = collect(6.0:1.0:20.0)
day       = 172.0

# Accumulate hourly W/m² into a daily Wh/m² total via the extractor:
# (run once per pixel across all hours, trapezoidal integration)
result = solar_radiation_grid(terrain_grids, solar_model, [day], hours_vec;
    extract = r -> ustrip.(u"W/m^2", r.global_terrain))
# result is (ny, nx, nhours) — integrate externally, or use the built-in:
result_daily = solar_radiation_grid(terrain_grids, solar_model, [day], hours_vec)
```
"""
function solar_radiation_grid(
    terrain_grids,
    solar_model,
    days,
    hours;
    albedo                         = 0.2,
    extract::F                     = r -> ustrip(u"W/m^2", r.global_terrain[1]),
    verbose::Bool                  = true,
) where {F}
    (; elevation_m, slope_deg, aspect_deg, latitude_deg, longitude_deg,
       pressure_r, horizons_u, data_is_xy) = terrain_grids

    ny, nx = size(elevation_m, data_is_xy ? 2 : 1),
             size(elevation_m, data_is_xy ? 1 : 2)

    alb_grid = albedo isa AbstractMatrix ? albedo : fill(Float64(albedo), ny, nx)

    days_f  = collect(Float64, days)
    hours_f = collect(Float64, hours)

    # Probe one valid pixel to learn the output shape from extract
    sample_result = nothing
    for I in CartesianIndices((ny, nx))
        i, j   = I[1], I[2]
        ri, rj = data_is_xy ? (j, i) : (i, j)
        lat = latitude_deg[ri, rj];  lon = longitude_deg[ri, rj]
        elev = elevation_m[ri, rj];  slp = slope_deg[ri, rj]
        asp  = aspect_deg[ri, rj];   pres = pressure_r[ri, rj]
        any(ismissing.([lat, lon, elev, slp, asp, pres])) && continue
        st = SolarTerrain(;
            latitude = lat, longitude = lon, elevation = elev,
            slope = slp, aspect = asp, albedo = alb_grid[i, j],
            atmospheric_pressure = pres,
            horizon_angles = @view(horizons_u[i, j, :]),
        )
        raw = solar_radiation(solar_model; solar_terrain = st,
                              days = days_f, hours = hours_f)
        sample_result = extract(raw)
        break
    end
    isnothing(sample_result) && error("No valid pixels found in terrain grid.")

    if sample_result isa AbstractVector
        n  = length(sample_result)
        ET = eltype(sample_result)
        out = Array{Union{ET, Missing}}(missing, ny, nx, n)
    else
        ET = typeof(sample_result)
        out = Array{Union{ET, Missing}}(missing, ny, nx)
    end

    n_total = ny * nx
    n_done  = Threads.Atomic{Int}(0)

    Threads.@threads :static for I in CartesianIndices((ny, nx))
        i, j   = I[1], I[2]
        ri, rj = data_is_xy ? (j, i) : (i, j)

        lat  = latitude_deg[ri, rj];  lon  = longitude_deg[ri, rj]
        elev = elevation_m[ri, rj];   slp  = slope_deg[ri, rj]
        asp  = aspect_deg[ri, rj];    pres = pressure_r[ri, rj]
        alb  = alb_grid[i, j]

        if any(ismissing.([lat, lon, elev, slp, asp, pres])) || ismissing(alb)
            continue
        end

        st = SolarTerrain(;
            latitude             = lat,
            longitude            = lon,
            elevation            = elev,
            slope                = slp,
            aspect               = asp,
            albedo               = alb,
            atmospheric_pressure = pres,
            horizon_angles       = @view(horizons_u[i, j, :]),
        )

        raw = solar_radiation(solar_model; solar_terrain = st,
                              days = days_f, hours = hours_f)
        val = extract(raw)

        if val isa AbstractVector
            out[i, j, :] .= val
        else
            out[i, j] = val
        end

        d = Threads.atomic_add!(n_done, 1)
        if verbose && Threads.threadid() == 1 && (d % max(1, nx) == 0 || d == n_total - 1)
            pct = round(100 * (d + 1) / n_total; digits = 1)
            print("  $(d+1)/$n_total pixels ($pct%)   \r")
        end
    end

    verbose && println()
    return out
end
