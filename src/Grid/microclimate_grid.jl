# Per-pixel gridded microclimate simulation.
#
# Provides `simulate_microclimate_grid`, which runs `simulate_microclimate`
# for every pixel in a terrain grid, with pluggable warm-start strategies
# via `GridInitStrategy` subtypes.

# ============================================================================
# Initialisation strategies
# ============================================================================

"""
    GridInitStrategy

Abstract supertype for strategies that supply the initial soil-temperature
profile for each pixel in a gridded [`simulate_microclimate_grid`](@ref) run.

| Subtype | Description |
|---------|-------------|
| [`RowWarmStart`](@ref)     | Converged midnight profile from the pixel above |
| [`ERA5SoilInit`](@ref)     | Pre-fetched ERA5 (or similar) soil temperature field |
| [`SolarMatchInit`](@ref)   | Nearest already-simulated pixel by daily solar load |
| [`LatitudeMatchInit`](@ref)| Nearest-latitude already-simulated row |
"""
abstract type GridInitStrategy end

"""
    RowWarmStart()

Use the converged midnight soil-temperature profile of the pixel directly above
as the initial condition for the current pixel.

Falls back to an elevation-based linear estimate for the first row or when the
pixel above has no data.  Works best for grids processed row-by-row where
adjacent pixels share similar thermal regimes.
"""
struct RowWarmStart <: GridInitStrategy end

"""
    ERA5SoilInit(soil_temperature)

Initialise each pixel from a gridded soil temperature field (e.g. ERA5 `stl1`).

`soil_temperature` must be an `(ny, nx, ndepths)` array of temperature values
with Kelvin units.  Pixels where any depth is `NaN` fall back to the
elevation-based linear estimate.
"""
struct ERA5SoilInit{A <: AbstractArray} <: GridInitStrategy
    soil_temperature::A   # (ny, nx, ndepths), Kelvin-unit values
end

"""
    SolarMatchInit(daily_solar)

Initialise each pixel from the already-simulated pixel with the most similar
cumulative daily solar radiation load.

`daily_solar` is a `(ny, nx)` matrix of pre-computed daily solar totals (any
consistent unit, e.g. Wh/m²) as returned by the solar radiation grid pipeline.

Useful for fine-scale terrain grids where a pixel on a south-facing slope may
be thermally closer to a flat pixel at lower elevation than to the pixel
directly above it.  Falls back to [`RowWarmStart`](@ref) until at least one
pixel has been simulated.
"""
struct SolarMatchInit{A <: AbstractMatrix} <: GridInitStrategy
    daily_solar::A   # (ny, nx) pre-computed daily solar totals
end

"""
    LatitudeMatchInit()

Initialise each pixel from the nearest-latitude already-simulated pixel in the
same column, falling back to any completed pixel in the nearest row.

Useful for coarse-resolution grids processed in arbitrary order where the
dominant temperature gradient is meridional.
"""
struct LatitudeMatchInit <: GridInitStrategy end

# ============================================================================
# Grid simulation state
# ============================================================================

"""
    GridSimState

Mutable state updated after each pixel completes so that subsequent pixels can
warm-start from it.

Fields:
- `completed_profiles`: `(ny, nx)` matrix of converged midnight soil-temperature
  profiles (`Vector{<:Quantity}` or `nothing` for not-yet-simulated pixels).
- `n_done`: thread-safe counter of completed pixels.
"""
mutable struct GridSimState
    completed_profiles :: Matrix{Any}
    n_done             :: Threads.Atomic{Int}
end

GridSimState(ny::Int, nx::Int) = GridSimState(
    fill!(Matrix{Any}(undef, ny, nx), nothing),
    Threads.Atomic{Int}(0),
)

# ============================================================================
# Initial-profile helpers
# ============================================================================

# Linear estimate: upper nodes at Tmean, second-deepest at the midpoint,
# deepest node at Tdeep.
function _elevation_T_init(ndepths, Tmean, Tdeep)
    T = Vector{typeof(1.0u"K")}(undef, ndepths)
    n = max(1, ndepths - 2)
    T[1:n]       .= Tmean
    ndepths >= 2 && (T[ndepths - 1] = (Tmean + Tdeep) / 2)
    T[ndepths]    = Tdeep
    return T
end

function _get_T_init(::RowWarmStart, state::GridSimState,
                     i, j, ndepths, Tmean, Tdeep)
    i > 1 && !isnothing(state.completed_profiles[i - 1, j]) &&
        return state.completed_profiles[i - 1, j]
    return _elevation_T_init(ndepths, Tmean, Tdeep)
end

function _get_T_init(init::ERA5SoilInit, _state::GridSimState,
                     i, j, ndepths, Tmean, Tdeep)
    profile = [init.soil_temperature[i, j, d] for d in 1:ndepths]
    any(x -> isnan(ustrip(u"K", x)), profile) &&
        return _elevation_T_init(ndepths, Tmean, Tdeep)
    return profile
end

function _get_T_init(init::SolarMatchInit, state::GridSimState,
                     i, j, ndepths, Tmean, Tdeep)
    solar_ij = init.daily_solar[i, j]
    if !isnan(solar_ij)
        best_profile = nothing
        best_dist    = Inf
        for ii in axes(state.completed_profiles, 1),
            jj in axes(state.completed_profiles, 2)
            p = state.completed_profiles[ii, jj]
            isnothing(p) && continue
            # cumulative_solar was stored by the grid loop below
            # Access via the init struct's daily_solar as a proxy key:
            # compare against the source pixel's pre-computed solar
            d = abs(init.daily_solar[ii, jj] - solar_ij)
            d < best_dist && (best_profile = p; best_dist = d)
        end
        !isnothing(best_profile) && return best_profile
    end
    # Fallback: row warm-start
    i > 1 && !isnothing(state.completed_profiles[i - 1, j]) &&
        return state.completed_profiles[i - 1, j]
    return _elevation_T_init(ndepths, Tmean, Tdeep)
end

function _get_T_init(::LatitudeMatchInit, state::GridSimState,
                     i, j, ndepths, Tmean, Tdeep)
    # Same column, nearest completed row above
    for row in (i - 1):-1:1
        p = state.completed_profiles[row, j]
        !isnothing(p) && return p
    end
    # Any completed pixel in the nearest row above
    for row in (i - 1):-1:1
        for jj in axes(state.completed_profiles, 2)
            p = state.completed_profiles[row, jj]
            !isnothing(p) && return p
        end
    end
    return _elevation_T_init(ndepths, Tmean, Tdeep)
end

# ============================================================================
# Grid simulation
# ============================================================================

"""
    simulate_microclimate_grid(
        terrain_grids, wp_grid, soil_thermal_model, solar_model;
        init_strategy, snapshot_steps, depths, heights,
        depth_indices, height_indices, albedo,
        roughness_height, karman_constant, dyer_constant, viewfactor,
        iterate_day, vapour_pressure_equation, verbose
    ) -> (; soil_temperature, air_temperature)

Run [`simulate_microclimate`](@ref) for every pixel in a terrain grid.

Rows are processed sequentially; columns within each row are parallelised with
`Threads.@threads`.  Start Julia with `--threads auto` for best performance.

# Positional arguments
- `terrain_grids`      : NamedTuple from [`compute_terrain_grids`](@ref)
- `wp_grid`            : `(ny, nx)` `Matrix` of per-pixel weather structs
  (lapse-corrected using e.g. [`lapse_correct_weather`](@ref)); `nothing`
  entries are skipped.
- `soil_thermal_model` : shared [`CampbelldeVriesSoilThermal`](@ref) instance
- `solar_model`        : [`SolarProblem`](@ref) instance

# Keyword arguments
- `init_strategy`  : [`GridInitStrategy`](@ref) controlling how each pixel's
  initial soil-temperature profile is chosen (default [`RowWarmStart`](@ref))
- `snapshot_steps` : 1-based time-step indices to store (default `1:24`)
- `depths`         : soil node depths, passed to `simulate_microclimate`
- `heights`        : air-temperature heights, passed to `simulate_microclimate`
- `depth_indices`  : `Int` or `AbstractVector{Int}` — subset of depth indices
  to store; `nothing` keeps all depths
- `height_indices` : `Int` or `AbstractVector{Int}` — subset of height indices
  to store; `nothing` keeps all heights
- `albedo`         : scalar or `(ny, nx)` matrix (default `0.15`)
- `roughness_height`, `karman_constant`, `dyer_constant`, `viewfactor`
  : MicroTerrain parameters (scalar, default values match NicheMapR)
- `iterate_day`    : maximum diurnal convergence iterations (default `5`)
- `vapour_pressure_equation` : passed to `simulate_microclimate`
- `verbose`        : print row-by-row progress (default `true`)

# Returns
NamedTuple with:
- `soil_temperature` : `Array{Float64,4}` `(ny, nx, nsteps, ndepths)` in °C
- `air_temperature`  : `Array{Float64,4}` `(ny, nx, nsteps, nheights)` in °C

Pixels with missing terrain data or `nothing` weather are left as `NaN`.
"""
function simulate_microclimate_grid(
    terrain_grids,
    wp_grid::AbstractMatrix,
    soil_thermal_model,
    solar_model;
    init_strategy::GridInitStrategy = RowWarmStart(),
    snapshot_steps                   = 1:24,
    depths                           = DEFAULT_DEPTHS,
    heights                          = [0.01, 2.0]u"m",
    depth_indices  :: Union{Nothing, Int, AbstractVector{Int}} = nothing,
    height_indices :: Union{Nothing, Int, AbstractVector{Int}} = nothing,
    albedo                           = 0.15,
    roughness_height                 = 0.004u"m",
    karman_constant                  = 0.4,
    dyer_constant                    = 16.0,
    viewfactor                       = 1.0,
    iterate_day    :: Int            = 5,
    vapour_pressure_equation         = GoffGratch(),
    verbose        :: Bool           = true,
)
    (; elevation_m, slope_deg, aspect_deg, latitude_deg, longitude_deg,
       pressure_r, horizons_u, data_is_xy) = terrain_grids

    ny, nx   = size(wp_grid)
    ndepths  = length(depths)
    nheights = length(heights)
    nsteps   = length(snapshot_steps)

    di = isnothing(depth_indices)  ? (1:ndepths)  :
         depth_indices  isa Int    ? [depth_indices]  : depth_indices
    hi = isnothing(height_indices) ? (1:nheights) :
         height_indices isa Int    ? [height_indices] : height_indices

    alb_grid = albedo isa AbstractMatrix ? albedo : fill(Float64(albedo), ny, nx)

    state   = GridSimState(ny, nx)
    n_total = ny * nx

    soil_out = fill(NaN, ny, nx, nsteps, length(di))
    air_out  = fill(NaN, ny, nx, nsteps, length(hi))

    verbose && println("Running per-pixel microclimate ($ny × $nx pixels, " *
                       "$(Threads.nthreads()) thread(s))...")

    @time for i in 1:ny
        Threads.@threads :static for j in 1:nx
            ri, rj = data_is_xy ? (j, i) : (i, j)

            elev = elevation_m[ri, rj]
            lat  = latitude_deg[ri, rj]
            lon  = longitude_deg[ri, rj]
            slp  = slope_deg[ri,    rj]
            asp  = aspect_deg[ri,   rj]
            pres = pressure_r[ri,   rj]
            wp   = wp_grid[i, j]

            if ismissing(elev) || ismissing(lat) || ismissing(lon) ||
               ismissing(slp)  || ismissing(asp) || ismissing(pres) || isnothing(wp)
                continue
            end

            alb = alb_grid[i, j]
            (ismissing(alb) || isnan(alb)) && continue

            mm    = wp.environment_minmax
            ed    = wp.environment_daily
            Tmean = (mm.reference_temperature_min[1] + mm.reference_temperature_max[1]) / 2
            Tdeep = ed.deep_soil_temperature[1]

            T_init = _get_T_init(init_strategy, state, i, j, ndepths, Tmean, Tdeep)

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
            mt = MicroTerrain(;
                elevation        = elev,
                roughness_height,
                karman_constant,
                dyer_constant,
                viewfactor,
            )

            result = simulate_microclimate(
                st, mt, soil_thermal_model, wp;
                depths, heights, solar_model,
                initial_soil_temperature = T_init,
                vapour_pressure_equation,
                iterate_day,
            )

            # Save converged midnight profile for warm-starting the row below
            state.completed_profiles[i, j] = collect(result.soil_temperature[1, :])

            @inbounds for (k, s) in enumerate(snapshot_steps)
                for (odi, rdi) in enumerate(di)
                    soil_out[i, j, k, odi] =
                        ustrip(u"°C", result.soil_temperature[s, rdi])
                end
                for (ohi, rhi) in enumerate(hi)
                    air_out[i, j, k, ohi] =
                        ustrip(u"°C", result.profile[s].air_temperature[rhi])
                end
            end

            d = Threads.atomic_add!(state.n_done, 1)
            if verbose && Threads.threadid() == 1 &&
               (d % max(1, nx) == 0 || d == n_total - 1)
                pct = round(100 * (d + 1) / n_total; digits = 1)
                print("  row $i/$ny  ($(d+1)/$n_total, $pct%)   \r")
            end
        end
    end

    verbose && println("\nSimulation complete.")
    return (; soil_temperature = soil_out, air_temperature = air_out)
end
