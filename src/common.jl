# micro_common.jl
#
# Shared types and scaffolding for the two declarative microclimate-simulation
# entry points:
#
#   MicroRasterProblem (grid mode, src/micro_raster.jl)
#   MicroVectorProblem (points mode, src/micro_vector.jl)
#
# Both produce a `MicroMapCache` and run through the same `solve!` /
# `_solve_proto_pixel!` / `_solve_remaining!` loop — the only difference is
# the spatial layout of the forcings (a 2-D `(X, Y)` grid vs a 1-D
# `Dim{:point}` lookup carrying point coordinates).
#
# The mode-specific `init` methods load and spatially-extract their forcings
# (weather, terrain, albedo, roughness, canonical overrides) into matching
# `Raster` shapes, then hand off to the shared `_build_inputs_and_pool`
# factory which produces the per-pixel `build_inputs` closure and worker pool
# used by `solve!`.

# ---------------------------------------------------------------------------
# Output layer specs and per-pixel write
# ---------------------------------------------------------------------------

# A `LayerSpec{Name, Kind}` encodes both the result-field name and its storage
# kind in its type, so each unrolled iteration specialises on a concrete type.
# `Kind` is one of:
#   :soil    — result.<Name> is a (Ti, Depth) matrix
#   :profile — result.profile.<Name> is a (Ti, Height) matrix
#   :scalar  — result.<Name> is a (Ti,) vector
# Number of azimuth directions used for horizon-angle computation. A
# power of two makes the cardinal directions land on sample points.
const N_HORIZON_ANGLES = 32

struct LayerSpec{Name, Kind} end
LayerSpec(name::Symbol, kind::Symbol) = LayerSpec{name, kind}()

@inline _layer_name(::LayerSpec{N}) where N = N
@inline _layer_source(result, ::LayerSpec{N, :profile}) where N = getproperty(result.profile, N)
@inline _layer_source(result, ::LayerSpec{N}) where N = getproperty(result, N)

const _DEFAULT_OUTPUT_LAYERS = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:air_temperature, :profile),
    LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed, :profile),
    LayerSpec(:surface_water, :scalar),
    LayerSpec(:global_radiation, :scalar),
    LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:snow_depth, :scalar),
)

"""
    canonical_unit(name::Symbol)

Canonical `Unitful` unit each default output layer is reported in for
storage and downstream analysis. Used by `strip_to_canonical` when
preparing a stack for unit-free I/O (e.g. NetCDF).
"""
canonical_unit(name::Symbol)              = canonical_unit(Val(name))
canonical_unit(::Val{:soil_temperature})  = u"°C"
canonical_unit(::Val{:air_temperature})   = u"°C"
canonical_unit(::Val{:sky_temperature})   = u"°C"
canonical_unit(::Val{:snow_depth})        = u"cm"
canonical_unit(::Val{:wind_speed})        = u"m/s"
canonical_unit(::Val{:surface_water})     = u"kg/m^2"
canonical_unit(::Val{:global_radiation})  = u"W/m^2"
canonical_unit(::Val{:soil_moisture})     = u"m^3/m^3"
canonical_unit(::Val{:relative_humidity}) = u"percent"

"""
    strip_to_canonical(stack::RasterStack)

Convert every layer of `stack` to its `canonical_unit` and strip the
Unitful unit, returning a Float-valued `RasterStack` suitable for writing
to a format that does not carry Unitful metadata (e.g. NetCDF).
"""
function strip_to_canonical(stack)
    names = propertynames(stack)
    return RasterStack(NamedTuple{names}(map(names) do n
        ustrip.(canonical_unit(n), stack[n])
    end))
end

# ---------------------------------------------------------------------------
# Shared types
# ---------------------------------------------------------------------------

"""
    MicroMapModel(; micro_model, dem_source, weather_source,
                    landcover_source=nothing,
                    surface_albedo_source=nothing,
                    roughness_height_source=nothing,
                    output_layers=_DEFAULT_OUTPUT_LAYERS,
                    lapse_rate_model=EnvironmentalLapseRate(),
                    solar_only=false,
                    solar_output_layers=())

Declarative description of a microclimate run. Constant across spatial
extents and time ranges — pair with a `MicroRasterProblem` (grid) or
`MicroVectorProblem` (points) to actually run.

- `micro_model::MicroModel` — inner per-pixel physics model
- `dem_source` — DEM data-source type (e.g. `SRTM`)
- `weather_source` — weather data-source type (e.g. `TerraClimate{Historical}`)
- `landcover_source` — optional land-cover data-source type
  (e.g. `EarthEnv{LandCover}`)
- `surface_albedo_source`, `roughness_height_source` — accepted forms:
    * `nothing` (use `landcover_source`'s default class→value table)
    * a class→value `NamedTuple` (overrides the dataset default)
    * a scalar (broadcast)
    * a `Raster` (resampled/extracted to the run template)
    * a `Type{<:RasterDataSource}` (load a property raster directly)
- `output_layers` — tuple of `LayerSpec`s controlling which output rasters
  are materialised
- `soil_moisture_source` — optional soil moisture data-source type
  (e.g. `CPCSoil`). Loaded automatically at `init` time and used as
  time-varying prescribed soil moisture. Overridden by `data.soil_moisture`
  if both are supplied.
- `lapse_rate_model::LapseRate` — atmospheric lapse-rate model for
  elevation-correcting weather data
- `solar_only::Bool` — when `true`, skip the microclimate ODE and return
  only solar radiation output. Defaults to four broadband layers when
  `solar_output_layers` is empty.
- `solar_output_layers` — tuple of `SolarOutputLayer`s controlling which
  solar radiation rasters are written alongside (or instead of) microclimate
  output. Each layer may request broadband or a waveband integral. Predefined
  constants: `SOLAR_BROADBAND`, `SOLAR_PAR`, `SOLAR_UVB`, `SOLAR_NIR`.
"""
@kwdef struct MicroMapModel{MM,DS,WS,LCS,SAS,RHS,SMS,OL,LRT,SOL}
    micro_model::MM
    dem_source::DS
    weather_source::WS
    landcover_source::LCS = nothing
    surface_albedo_source::SAS = nothing
    roughness_height_source::RHS = nothing
    soil_moisture_source::SMS = nothing
    output_layers::OL = _DEFAULT_OUTPUT_LAYERS
    lapse_rate_model::LRT = EnvironmentalLapseRate()
    compute_terrain::Bool = true
    solar_only::Bool = false
    cloud_correct_solar::Bool = false
    solar_output_layers::SOL = ()
end

# Per-run workspace built by `init(problem)`. Holds the spatially-extracted
# forcings, pre-resolved surface grids, and a Channel pool of per-worker
# inner `MicroCache`s + scratch buffers. The output `RasterStack` is
# allocated fresh on every `solve!` call — it does not live on the cache.
#
# Used by both `MicroRasterProblem` and `MicroVectorProblem`: only the shape of
# the spatial dim differs (`(X, Y)` vs `(Dim{:point},)`); every solve-time
# code path indexes via `I::Tuple` of dim wrappers from `DimIndices`, so the
# loop body is mode-agnostic.
mutable struct MicroMapCache{P,W,T,A,R,CO,M,POOL,SC,CC}
    problem::P                       # MicroRasterProblem/MicroVectorProblem
    weather::W                       # RasterStack
    terrain::T                       # RasterStack
    albedo_grid::A                   # Raster
    roughness_grid::R                # Raster
    canonical_overrides::CO          # NamedTuple of pre-extracted Rasters
    mask::M                          # nothing, or a 2-D Bool Raster (X, Y)
    cache_pool::POOL                 # Channel{(; micro::MicroCache, scratch)}
    init_inputs::SC                  # Resolved initial conditions + build_inputs closure
    cloud_constants::CC              # Shared SolarRadiation constants
end

"""
    terrain(cache::MicroMapCache) -> RasterStack

Return the terrain `RasterStack` computed during `init`. Layers:
`elevation`, `slope`, `aspect`, `latitude`, `longitude`,
`atmospheric_pressure`, `horizon_angles`.

Pass this as `data = (; terrain = terrain(cache))` on a subsequent
`MicroRasterProblem` or `MicroVectorProblem` to reuse it and skip the
DEM download and horizon-angle computation.
"""
terrain(cache::MicroMapCache) = cache.terrain

# ---------------------------------------------------------------------------
# Init helpers — shared by both modes
# ---------------------------------------------------------------------------

# Keys in `problem.data` that target collection sources or named single
# layers; anything else is interpreted as a canonical weather variable.
const _SPECIAL_DATA_KEYS = (
    :weather, :landcover, :dem, :terrain, :surface_albedo, :roughness_height,
)

@inline _is_special_key(k::Symbol) = k in _SPECIAL_DATA_KEYS

# Filter `data` down to the canonical-weather-variable overrides
# (anything not in `_SPECIAL_DATA_KEYS`).
function _canonical_data(data::NamedTuple)
    ks = filter(k -> !_is_special_key(k), keys(data))
    return NamedTuple{ks}(map(k -> data[k], ks))
end

# Resolve initial conditions. Every key defaults to `nothing`, meaning
# "fall back to the weather source / data override / inner-solver default"
# — never fabricated. User-supplied keys override; missing keys are filled
# from `_DEFAULT_INIT` so the build_inputs closure always sees a complete
# NamedTuple.
const _DEFAULT_INIT = (
    soil_temperature = nothing,
    soil_moisture = nothing,
    snow_depth = nothing,
    snow_temperature = nothing,
    snow_density = nothing,
)
_resolve_init(it::NamedTuple, _) = merge(_DEFAULT_INIT, it)
_resolve_init(::Nothing, _) = _DEFAULT_INIT

# Does the run actually have a `:soil_moisture` value to read? True iff
# the weather source declares it as a canonical variable or the user
# supplied `data.soil_moisture`.
function _has_canonical_input(name::Symbol, weather_source, data::NamedTuple)
    haskey(data, name) && return true
    weather_source === nothing && return false
    for var in variables(weather_source)
        canonical_name(var) === name && return true
    end
    return false
end

# Per-pixel `initial_soil_moisture` resolution.
#  - User-supplied vector wins.
#  - Otherwise read the first timestep of `buffers.soil_moisture` and
#    broadcast across depths (TerraClimate-style rooting-zone average).
#  - If neither is available, refuse to fabricate.
_initial_soil_moisture(user::AbstractVector, _, _, _) = user
function _initial_soil_moisture(::Nothing, buffer, available::Bool, depths)
    available || error(
        "No soil_moisture available. Either: pass `init = (; soil_moisture = vec)`, " *
        "use a weather source that provides it as a canonical variable (e.g. TerraClimate), " *
        "or supply `data.soil_moisture` on the problem."
    )
    return fill(buffer[1], length(depths))
end

# Source-loading helpers (override beats loader). Mode-agnostic: both grid
# and points init build their own `area::Extent` (rectangular for grid;
# bounding box of points for points mode) and call these.
_resolve_dem(data::NamedTuple, source, area) =
    haskey(data, :dem) ? data.dem : _load_dem(source, area)

_resolve_weather(data::NamedTuple, source, area, years) =
    haskey(data, :weather) ? data.weather : _load_weather(source, area, years)

# Normalise `area` to an Extent for loaders; geometry is kept separately for the mask.
_to_extent(area::Extent) = area
_to_extent(area) = GeoInterface.extent(area)


# ---------------------------------------------------------------------------
# Source-label helper — shared by raster and vector init @info messages
# ---------------------------------------------------------------------------

function _soil_moisture_label(strategy, soil_moisture_available, init_inputs)
    if nameof(typeof(strategy)) === :PrescribedSoilMoisture
        hasproperty(strategy, :precomputed_soil_moisture) &&
            !isnothing(strategy.precomputed_soil_moisture) && return "prescribed (time-varying, precomputed)"
        soil_moisture_available                            && return "prescribed (time-varying, from weather source)"
        !isnothing(init_inputs.soil_moisture)             && return "prescribed (fixed, user-supplied)"
        return "prescribed (fixed, default)"
    end
    return string(nameof(typeof(strategy)))
end

_source_label(::Nothing)                                      = "from landcover"
_source_label(x::Number)                                      = "fixed: $x"
_source_label(x::Quantity)                                    = "fixed: $x"
_source_label(::Raster)                                       = "user Raster"
_source_label(s::Type{<:RasterDataSources.RasterDataSource})  = string(s)

# ---------------------------------------------------------------------------
# Native surface-property resolution
# ---------------------------------------------------------------------------
#
# `_resolve_surface_native` returns a *native-resolution* representation of a
# surface property (albedo / roughness) — either a scalar (for the "broadcast
# everywhere" case) or a 2-D `(X, Y)` `Raster`. Mode-specific code then runs
# this through the appropriate spatial extractor (`_resample(...; to=template)`
# in grid mode, per-point `Near`-lookup in points mode).
#
# The forms accepted by `MicroMapModel.surface_albedo_source` /
# `.roughness_height_source` are identical to the legacy single-step path:
#   - `nothing` → dispatch landcover-source default class→value table
#   - `NamedTuple` of weights → load landcover, weight at native res
#   - scalar → return the scalar; mode wraps it in a `Fill`-Raster
#   - `Raster` → return as-is
#   - `Type{<:RasterDataSource}` → load via `load_surface_property`

# Scalar / unrecognised value: return as-is. Mode extractor builds the Fill.
_resolve_surface_native(value, _landcover_source, _area, _default_fn) = value
# User-supplied `Raster`: return as-is. Mode extractor resamples/extracts.
_resolve_surface_native(r::Raster, _landcover_source, _area, _default_fn) = r
# `nothing` → dispatch the default class→value NamedTuple for the
# landcover source, then re-enter via the NamedTuple branch.
_resolve_surface_native(::Nothing, landcover_source, area, default_fn) =
    _resolve_surface_native(default_fn(landcover_source), landcover_source, area, default_fn)
_resolve_surface_native(::Nothing, ::Nothing, _area, default_fn) =
    error("`$(default_fn)`: no value supplied. Pass a `landcover_source` " *
          "(e.g. `EarthEnv{LandCover}`, `MODIS{MCD12Q1}`), a `Type{<:RasterDataSource}` " *
          "for the property itself, or a scalar / Raster / NamedTuple.")
# Class→value NamedTuple: load landcover via `landcover_source` and weight
# at native resolution. Returned as a `Raster` carrying `commondims(landcover, (X, Y))`.
function _resolve_surface_native(
    weights::NamedTuple,
    landcover_source::Type{<:RasterDataSources.RasterDataSource},
    area, _default_fn,
)
    landcover = load_landcover(landcover_source, area)
    weighted_native = landcover_weighted(landcover, weights, landcover_source)
    spatial_dims = commondims(landcover, (X(), Y()))
    return _wrap_native_raster(weighted_native, spatial_dims)
end
# `Type{<:RasterDataSource}` as the property *source* itself (not landcover):
# load a 2-D Raster from that dataset; ignores `landcover_source`. Each
# dataset that supports this implements `load_surface_property(source, area)`.
_resolve_surface_native(
    property_source::Type{<:RasterDataSources.RasterDataSource},
    _landcover_source, area, _default_fn,
) = load_surface_property(property_source, area)

@inline _wrap_native_raster(values::Raster, _spatial_dims) = values
@inline _wrap_native_raster(values, spatial_dims) = Raster(values, spatial_dims)

# ---------------------------------------------------------------------------
# Build the per-pixel inputs closure + worker pool
# ---------------------------------------------------------------------------

# Build the cloud-derivation constants. Immutable, build-once-per-run.
function _build_cloud_constants()
    return (;
        solar_model = SolarProblem(),
        hours = collect(0.0:1.0:23.0),
        flat_terrain_template = SolarTerrain(;
            elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
            horizon_angles = Fill(0.0u"°", N_HORIZON_ANGLES),
            albedo = 0.15,
            atmospheric_pressure = atmospheric_pressure(0.0u"m"),
            latitude = 0.0u"°", longitude = 0.0u"°",
        ),
    )
end

# Build the `build_inputs(scratch, I)` closure and the per-worker cache pool
# shared by both grid and points init. All spatial forcings must already be
# laid out so that `forcing[I...]` works for every `I` from
# `DimIndices(terrain.elevation)` (i.e. `(X(i), Y(j))` in grid mode and
# `(Dim{:point}(p),)` in points mode).
function _build_inputs_and_pool(;
    model, weather_source, weather, terrain, mask,
    albedo_grid, roughness_grid, canonical_overrides,
    init_inputs, soil_moisture_available, years, days, cloud_constants,
    soil_profile, target_timestep::Timestep = Hourly(),
)
    (; micro_model, lapse_rate_model) = model
    vapour_pressure_method = micro_model.vapour_pressure_equation

    calendar = weather_calendar(weather_source)
    # Must match `days`, not the full `years` span, or sub-yearly runs
    # size weather buffers for a whole year while solar geometry only
    # covers the request.
    solar_ndays = length(days)
    time_mode = _time_mode(calendar)
    nsteps = length(cloud_constants.hours) * solar_ndays
    nmax = cloud_constants.solar_model.wavelength_count
    allocate_scratch() = (;
        weather = allocate_weather_buffers(weather_source, target_timestep, years; days_of_year = days),
        solar = (;
            out = allocate_output_arrays(nsteps, solar_ndays, nmax),
            buffers = allocate_buffers(nmax, cloud_constants.solar_model.diffuse_model),
        ),
        cloud_constants,
    )

    npixels = length(terrain.elevation)
    nworkers = min(Threads.nthreads(), npixels)

    if model.solar_only
        @info "model: threads: $(Threads.nthreads())"
        # No ODE needed — pool contains scratch-only workers (no micro field).
        # build_inputs is never called in _solve_solar_only!.
        _build_inputs_noop = (scratch, I) -> nothing
        proto = (; scratch = allocate_scratch())
        cache_pool = Channel{typeof(proto)}(nworkers)
        put!(cache_pool, proto)
        for _ in 2:nworkers
            put!(cache_pool, (; scratch = allocate_scratch()))
        end
        return (_build_inputs_noop, cache_pool)
    end

    wind_tgt = round(ustrip(u"m", maximum(micro_model.heights)); digits = 2)
    @info "model: snow:            $(nameof(typeof(micro_model.snow_model)))"
    sm_src = model.soil_moisture_source !== nothing ? " ($(model.soil_moisture_source))" : ""
    @info "model: soil moisture:   $(_soil_moisture_label(micro_model.config.soil_moisture_strategy, soil_moisture_available, init_inputs))$(sm_src)"
    @info "model: lapse rate:      $(nameof(typeof(lapse_rate_model)))"
    @info "model: wind:            reference height $(wind_tgt) m, power law-corrected from 10 m (source)"
    @info "model: threads:         $(Threads.nthreads())"

    function build_inputs(scratch, I::Tuple)
        horizon_angles = terrain.horizon_angles[I...]
        site = Site(;
            elevation = terrain.elevation[I...],
            slope = terrain.slope[I...],
            aspect = terrain.aspect[I...],
            latitude = terrain.latitude[I...],
            longitude = terrain.longitude[I...],
            atmospheric_pressure = terrain.atmospheric_pressure[I...],
            horizon_angles,
            sky_view_fraction = _sky_view_from_horizon(horizon_angles),
            albedo = albedo_grid[I...],
            roughness_height = roughness_grid[I...],
        )
        # TODO: why is this hacked in - there are two elevations
        ge = weather_grid_elevation(weather_source, weather, I)
        env = assemble_weather!(scratch, weather, weather_source, site, I;
            vapour_pressure_method, lapse_rate_model, canonical_overrides,
            grid_elevation = isnothing(ge) ? site.elevation : ge,
            wind_reference_height = maximum(micro_model.heights),
            target_timestep)
        initial_soil_moisture = _initial_soil_moisture(
            init_inputs.soil_moisture, scratch.weather.soil_moisture,
            soil_moisture_available, micro_model.depths,
        )
        # Snow init: fall back to MicroInputs' own defaults when the user
        # didn't supply a value. `initial_snow_density = nothing` is the
        # MicroInputs default ("use the snow model's snow_density") — pass
        # through unchanged.
        MicroInputs(; site, soil_profile,
            env.environment_minmax, env.environment_daily, env.environment_hourly,
            initial_soil_temperature = init_inputs.soil_temperature,
            initial_soil_moisture,
            initial_snow_depth = something(init_inputs.snow_depth, 0.0u"cm"),
            initial_snow_temperature = something(init_inputs.snow_temperature, u"K"(0.0u"°C")),
            initial_snow_density = init_inputs.snow_density,
        )
    end

    ci = findfirst(mask)
    isnothing(ci) &&
        error("All pixels are masked or have missing weather data (ocean?).")
    first_I = DimIndices(terrain.elevation)[ci]
    npixels = length(terrain.elevation)
    build_cache() = let scratch = allocate_scratch()
        (micro = CommonSolve.init(MicroProblem(micro_model, build_inputs(scratch, first_I); days, time_mode)),
         scratch)
    end

    proto = build_cache()
    cache_pool = Channel{typeof(proto)}(nworkers)
    put!(cache_pool, proto)
    for _ in 2:nworkers
        put!(cache_pool, build_cache())
    end

    return (build_inputs, cache_pool)
end

# ---------------------------------------------------------------------------
# CommonSolve.solve! / solve
# ---------------------------------------------------------------------------

"""
    CommonSolve.solve!(cache::MicroMapCache) -> RasterStack
    CommonSolve.solve!(output::RasterStack, cache::MicroMapCache) -> RasterStack

Run the per-pixel microclimate solve across every spatial index of
`cache.terrain.elevation`. In grid mode that's every `(X, Y)` pixel; in
points mode it's every `Dim{:point}` entry. Work is dispatched over the
cache pool with one `Threads.@spawn` worker per pool slot.

The one-arg form allocates a fresh `RasterStack` to write into — sized
from the first pixel's result so the output's element types are concrete.
The two-arg form writes into the supplied `output` and returns it, for
callers that want to reuse storage across runs.
"""
function CommonSolve.solve!(cache::MicroMapCache)
    model = cache.problem.model
    solar_pairs = cache.init_inputs.solar_pairs

    if model.solar_only
        @info "solve: solar-only mode ($(length(solar_pairs)) layer(s))..."
        solar_output = _allocate_solar_output(
            cache.terrain, solar_pairs, cache.init_inputs.anchor_dates, cache.mask)
        return _solve_solar_only!(solar_output, cache, solar_pairs)
    end

    @info "solve: solving first pixel..."
    proto, first_I, first_result = _solve_proto_pixel!(cache)
    layers = model.output_layers
    try
        output = _allocate_output(model.micro_model,
            cache.terrain, first_result, layers,
            cache.init_inputs.anchor_dates, cache.mask)
        _write_output!(output, first_result, layers, first_I)
        solar_output = if !isempty(solar_pairs)
            so = _allocate_solar_output(
                cache.terrain, solar_pairs, cache.init_inputs.anchor_dates, cache.mask)
            _compute_solar_for_pixel!(proto.scratch, cache.terrain, cache.albedo_grid, first_I)
            _write_solar_output!(so, proto.scratch.solar.out, solar_pairs,
                cache.cloud_constants.solar_model.wavelengths, first_I)
            so
        else
            nothing
        end
        @info "solve: starting main loop..."
        return _solve_remaining!(output, solar_output, cache, proto, first_I)
    catch
        put!(cache.cache_pool, proto)
        rethrow()
    end
end

function CommonSolve.solve!(output::RasterStack, cache::MicroMapCache)
    model = cache.problem.model
    solar_pairs = cache.init_inputs.solar_pairs
    proto, first_I, first_result = _solve_proto_pixel!(cache)
    layers = model.output_layers
    try
        _write_output!(output, first_result, layers, first_I)
        solar_output = if !isempty(solar_pairs)
            so = _allocate_solar_output(
                cache.terrain, solar_pairs, cache.init_inputs.anchor_dates, cache.mask)
            _compute_solar_for_pixel!(proto.scratch, cache.terrain, cache.albedo_grid, first_I)
            _write_solar_output!(so, proto.scratch.solar.out, solar_pairs,
                cache.cloud_constants.solar_model.wavelengths, first_I)
            so
        else
            nothing
        end
        return _solve_remaining!(output, solar_output, cache, proto, first_I)
    catch
        put!(cache.cache_pool, proto)
        rethrow()
    end
end

# Pull worker #1 from the pool and solve the first mask-active pixel that
# actually converges. This pixel sizes the output arrays, so we must find one
# that completes. Coastal pixels with post-resample NaN in derived quantities
# can still fail the ODE; we skip those silently and try the next candidate.
function _solve_proto_pixel!(cache)
    mask  = cache.mask
    proto = take!(cache.cache_pool)
    for I in DimIndices(cache.terrain.elevation)
        mask[I...] || continue
        reinit!(proto.micro, cache.init_inputs.build_inputs(proto.scratch, I))
        try
            solve!(proto.micro)
            return proto, I, proto.micro.output
        catch
        end
    end
    put!(cache.cache_pool, proto)
    error("No pixel solved successfully — all pixels are masked, ocean, or have invalid weather data.")
end

function _solve_remaining!(output, solar_output, cache, proto, first_I)
    cache_pool = cache.cache_pool
    layers = cache.problem.model.output_layers
    build_inputs = cache.init_inputs.build_inputs
    mask = cache.mask
    has_solar = solar_output !== nothing
    solar_pairs = cache.init_inputs.solar_pairs
    wavelengths = has_solar ? cache.cloud_constants.solar_model.wavelengths : nothing
    put!(cache_pool, proto)

    pixel_indices = DimIndices(cache.terrain.elevation)
    work = Channel{eltype(pixel_indices)}(max(length(pixel_indices) - 1, 1))
    nwork = 0
    @info "solve: building work channel ($(length(pixel_indices)) pixels)..."
    for I in pixel_indices
        I == first_I && continue
        mask[I...] || continue
        put!(work, I)
        nwork += 1
    end
    close(work)
    npixels_active = nwork + 1  # proto pixel already solved

    # Show progress only for non-trivial runs (> 10 active pixels).
    # The time check runs inside the worker loop (not a separate task) so it
    # fires even when all threads are saturated with ODE solving.
    show_progress = npixels_active > 10
    done = Threads.Atomic{Int}(1)
    t_start = time()
    last_report_ms = Threads.Atomic{Int64}(round(Int64, (t_start - 10.1) * 1000))

    nworkers = cache_pool.sz_max
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            c = take!(cache_pool)
            try
                for I in work
                    reinit!(c.micro, build_inputs(c.scratch, I))
                    try
                        solve!(c.micro)
                    catch
                        continue  # leave output as NaN for failed pixels
                    end
                    _write_output!(output, c.micro.output, layers, I)
                    if has_solar
                        _compute_solar_for_pixel!(c.scratch, cache.terrain, cache.albedo_grid, I)
                        _write_solar_output!(solar_output, c.scratch.solar.out,
                            solar_pairs, wavelengths, I)
                    end
                    n = Threads.atomic_add!(done, 1) + 1
                    if show_progress
                        t_now = time()
                        t_now_ms = round(Int64, t_now * 1000)
                        prev_ms = last_report_ms[]
                        if t_now_ms - prev_ms >= 10_000 &&
                                Threads.atomic_cas!(last_report_ms, prev_ms, t_now_ms) == prev_ms
                            elapsed = t_now - t_start
                            pct   = round(Int, 100 * n / npixels_active)
                            eta_s = round(Int, elapsed / n * (npixels_active - n))
                            @info "Raster solve: $n / $npixels_active pixels ($pct%) — ETA $(eta_s)s"
                        end
                    end
                end
            finally
                put!(cache_pool, c)
            end
        end
    end
    solar_output === nothing && return output
    return RasterStack(merge(NamedTuple(output), NamedTuple(solar_output)))
end

# Populate scratch.solar.out for pixel I using actual terrain geometry.
# Called for every solar output pass — both solar_only and combined modes.
# For sources whose derivation chain includes :solar_geometry (e.g. NCEP),
# this overwrites the flat-terrain solar computed in assemble_weather! with
# the slope/aspect/horizon-corrected values needed for global_terrain output.
# For monthly/daily sources (e.g. CRUCL2, TerraClimate), this is the only
# call that populates scratch.solar.out, which their derivation chain skips.
function _compute_solar_for_pixel!(scratch, terrain, albedo_grid, I)
    horizon_angles = terrain.horizon_angles[I...]
    solar_terrain = SolarTerrain(;
        elevation            = terrain.elevation[I...],
        slope                = terrain.slope[I...],
        aspect               = terrain.aspect[I...],
        latitude             = terrain.latitude[I...],
        longitude            = terrain.longitude[I...],
        atmospheric_pressure = terrain.atmospheric_pressure[I...],
        horizon_angles,
        albedo               = albedo_grid[I...],
    )
    solar_radiation!(scratch.solar.out, scratch.solar.buffers,
        scratch.cloud_constants.solar_model;
        solar_terrain,
        days  = scratch.weather.days_of_year,
        hours = scratch.cloud_constants.hours)
    return nothing
end

# Return the native weather field name that maps to canonical :cloud_cover,
# or `nothing` when the source does not provide cloud cover.
# Return the Variable for :cloud_cover (carries native field name +
# transform), or nothing if the source does not provide cloud cover.
function _cloud_weather_variable(weather_source)
    for var in variables(weather_source)
        canonical_name(var) === :cloud_cover && return var
    end
    return nothing
end

# Fill `factors` (length nsteps) with per-step Ångström sunshine fractions.
# Each weather Ti step covers `nhours_per_step` consecutive solar output steps
# (24 for monthly/daily sources). The Variable transform + unit are
# applied to the raw native field value before clamping to [0, 1] — this
# mirrors _copy_weather_to_buffers! so the cloud fraction is consistent
# regardless of whether the native field is already in [0,1] or needs
# conversion (e.g. CRUCL2 stores sunshine % via transform (100-s)/100).
function _fill_cloud_factors!(factors, weather, cloud_var::Variable, I, nhours_per_step)
    cloud_layer = getproperty(weather, native_field(cloud_var))
    n_weather = length(lookup(cloud_layer, Ti))
    k = 1
    for d in 1:n_weather
        raw   = Float64(cloud_layer[I..., Ti(d)])
        cloud = clamp(cloud_var.transform(raw) * cloud_var.unit, 0.0, 1.0)
        sf    = sunshine_fraction(Angstrom(), cloud)
        for _ in 1:nhours_per_step
            factors[k] = sf
            k += 1
        end
    end
    return factors
end

function _solve_solar_only!(solar_output, cache, solar_pairs)
    cache_pool = cache.cache_pool
    terrain = cache.terrain
    albedo_grid = cache.albedo_grid
    mask = cache.mask
    wavelengths = cache.cloud_constants.solar_model.wavelengths
    model = cache.problem.model

    # Cloud correction setup — only when weather was loaded and source has cloud cover.
    cloud_var = _cloud_weather_variable(model.weather_source)
    cloud_correct = model.cloud_correct_solar &&
        cloud_var !== nothing && !isempty(cache.weather)
    nsteps = size(first(values(solar_output)), Ti)
    nhours_per_step = if cloud_correct
        n_weather = length(lookup(getproperty(cache.weather, native_field(cloud_var)), Ti))
        nsteps ÷ n_weather
    else
        0
    end

    pixel_indices = DimIndices(terrain.elevation)
    work = Channel{eltype(pixel_indices)}(max(length(pixel_indices), 1))
    nwork = 0
    for I in pixel_indices
        mask[I...] || continue
        put!(work, I)
        nwork += 1
    end
    close(work)

    show_progress = nwork > 10
    done = Threads.Atomic{Int}(0)
    t_start = time()
    last_report_ms = Threads.Atomic{Int64}(round(Int64, (t_start - 10.1) * 1000))

    nworkers = cache_pool.sz_max
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            c = take!(cache_pool)
            # Pre-allocate one cloud-factor buffer per worker (reused across pixels).
            cloud_factors = cloud_correct ? Vector{Float64}(undef, nsteps) : nothing
            try
                for I in work
                    if cloud_correct
                        _fill_cloud_factors!(cloud_factors, cache.weather,
                                            cloud_var, I, nhours_per_step)
                    end
                    _compute_solar_for_pixel!(c.scratch, terrain, albedo_grid, I)
                    _write_solar_output!(solar_output, c.scratch.solar.out,
                        solar_pairs, wavelengths, I; cloud_factors)
                    if show_progress
                        n = Threads.atomic_add!(done, 1) + 1
                        t_now = time()
                        t_now_ms = round(Int64, t_now * 1000)
                        prev_ms = last_report_ms[]
                        if t_now_ms - prev_ms >= 10_000 &&
                                Threads.atomic_cas!(last_report_ms, prev_ms, t_now_ms) == prev_ms
                            elapsed = t_now - t_start
                            pct   = round(Int, 100 * n / nwork)
                            eta_s = round(Int, elapsed / n * (nwork - n))
                            @info "Solar solve: $n / $nwork pixels ($pct%) — ETA $(eta_s)s"
                        end
                    end
                end
            finally
                put!(cache_pool, c)
            end
        end
    end
    return solar_output
end

# ---------------------------------------------------------------------------
# Output stack allocation and per-pixel writes
# ---------------------------------------------------------------------------

function _allocate_output(model::MicroModel, terrain, proto, layers::Tuple,
                          anchor_dates::AbstractVector{Date}, mask)
    spatial_dims = dims(terrain.elevation)
    ti = Ti(_ti_datetime_axis(anchor_dates, model.hours))
    # Dim lookups must be plain numbers, not Unitful — Unitful in dims
    # breaks NetCDF I/O and surprises selectors. Strip both to metres so
    # `depth` and `height` share a single linear axis.
    extra = (
        soil = (ti, Dim{:depth}(ustrip.(u"m", model.depths))),
        profile = (ti, Dim{:height}(ustrip.(u"m", model.heights))),
        scalar = (ti,),
    )
    return RasterStack(NamedTuple(map(layers) do spec
        _layer_name(spec) => _allocate_layer(proto, spec, spatial_dims, extra, mask)
    end))
end

# Build the Ti axis as actual DateTimes. `anchor_dates` has one entry per
# solver step (one per unique day for daily/hourly; one per selected month
# for monthly). `hours` are 0..23 floats. The result is one DateTime per
# (step, hour) pair. Using actual calendar dates (not doy arithmetic) keeps
# cross-year and leap-year boundaries correct.
function _ti_datetime_axis(anchor_dates::AbstractVector{Date}, hours::AbstractVector)
    out = Vector{DateTime}(undef, length(anchor_dates) * length(hours))
    k = 1
    for d in anchor_dates, h in hours
        out[k] = DateTime(d) + Hour(round(Int, h))
        k += 1
    end
    return out
end

function _allocate_layer(proto, spec::LayerSpec{<:Any, K}, spatial_dims, extra, mask) where K
    T = typeof(first(_layer_source(proto, spec)))
    ds = (spatial_dims..., extra[K]...)
    return _allocate_layer_array(T, ds, mask)
end

_allocate_layer_array(T, ds, ::Nothing) = Raster(fill(_nan_of(T), map(length, ds)...), ds)
_allocate_layer_array(T, ds, _mask) =
    Rasters.create(nothing, T, ds; missingval = missing, fill = missing)

_nan_of(::Type{T}) where T<:AbstractFloat  = T(NaN)
_nan_of(::Type{T}) where T<:Quantity       = NaN * oneunit(T)

# `I` is a tuple of dim selectors from `DimIndices`: `(X(i), Y(j))` in grid
# mode, `(Dim{:point}(p),)` in points mode. The spatial axes are addressed
# by name so the storage order of `rast` is irrelevant. The trailing
# `Ti`/extra-dim selectors target the time and depth/height axes we
# constructed in `_allocate_layer`. Dispatch is on the **source's** rank
# (Vector = scalar, Matrix = soil/profile) so the same methods serve grid
# (rast 3-D / 4-D) and points (rast 2-D / 3-D) mode.
@inline function _write_slice!(rast, src::AbstractVector, I::Tuple)
    @inbounds for t in eachindex(src)
        rast[I..., Ti(t)] = src[t]
    end
    return nothing
end
@inline function _write_slice!(rast, src::AbstractMatrix, I::Tuple)
    # The trailing dim is `Dim{:depth}` for soil layers and `Dim{:height}`
    # for profile layers — extract its base type from the raster itself so
    # this method handles both kinds.
    extra_dim = basetypeof(last(dims(rast)))
    @inbounds for d in axes(src, 2), t in axes(src, 1)
        rast[I..., Ti(t), extra_dim(d)] = src[t, d]
    end
    return nothing
end

@inline function _write_output!(output, result, layers::Tuple, I::Tuple)
    unrolled_map(layers) do spec
        _write_slice!(output[_layer_name(spec)], _layer_source(result, spec), I)
        nothing
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Show helpers (shared)
# ---------------------------------------------------------------------------

# Each forcing has an origin describing where its value comes from at solve
# time. `_forcing_origin` classifies the origin from the model alone;
# `_forcing_origins(problem)` overlays the problem's `data` overrides on top.
_forcing_origin(::Val{:dem}, model) = _source_origin(model.dem_source, "data.dem")
_forcing_origin(::Val{:weather}, model) = _source_origin(model.weather_source, "data.weather")
_forcing_origin(::Val{:landcover}, model) =
    model.landcover_source === nothing ?
        (label = "(none)", missing = false) :
        (label = string(model.landcover_source), missing = false)
_forcing_origin(::Val{:surface_albedo}, model) =
    _surface_property_origin(model.surface_albedo_source, model.landcover_source, :surface_albedo)
_forcing_origin(::Val{:roughness_height}, model) =
    _surface_property_origin(model.roughness_height_source, model.landcover_source, :roughness_height)

# Apply problem.data overrides over the model-level origins. An override
# for a special key replaces the model's label with "user override".
function _forcing_origins(problem)
    data = problem.data
    map(_FORCINGS) do forcing
        origin = _forcing_origin(Val(forcing), problem.model)
        if haskey(data, forcing)
            return (forcing, (label = "user override (was: $(origin.label))", missing = false))
        end
        return (forcing, origin)
    end
end

_source_origin(::Nothing, data_key::String) =
    (label = "<missing — pass `$data_key`>", missing = true)
_source_origin(source, _) = (label = string(source), missing = false)

_surface_property_origin(::Nothing, ::Nothing, key) =
    (label = "<missing — pass landcover_source or `data.$key`>", missing = true)
_surface_property_origin(::Nothing, lc, _) =
    (label = "from $lc defaults", missing = false)
_surface_property_origin(::NamedTuple, ::Nothing, key) =
    (label = "<missing landcover_source for class weights>", missing = true)
_surface_property_origin(::NamedTuple, lc, _) =
    (label = "user weights via $lc", missing = false)
_surface_property_origin(s::Type{<:RasterDataSources.RasterDataSource}, _, _) =
    (label = "from $s", missing = false)
_surface_property_origin(::Raster, _, _) =
    (label = "user Raster", missing = false)
_surface_property_origin(v, _, _) =
    (label = "constant $v", missing = false)

const _FORCINGS = (:dem, :weather, :landcover, :surface_albedo, :roughness_height)

function Base.show(io::IO, ::MIME"text/plain", model::MicroMapModel)
    println(io, "MicroMapModel")
    origins = map(forcing -> (forcing, _forcing_origin(Val(forcing), model)), _FORCINGS)
    _print_forcing_table(io, origins)
    println(io, "  lapse_rate_model:  ", model.lapse_rate_model)
    println(io, "  output_layers:     ",
            join((_layer_name(spec) for spec in model.output_layers), ", "))
    missing_forcings = [string(forcing) for (forcing, o) in origins if o.missing]
    if !isempty(missing_forcings)
        println(io)
        println(io, "missing — pass via `data` on the problem:")
        for forcing in missing_forcings
            println(io, "  ", forcing)
        end
    end
end

function _print_forcing_table(io::IO, origins)
    width = maximum(length(string(forcing)) for (forcing, _) in origins)
    for (forcing, origin) in origins
        marker = origin.missing ? "⚠ " : "  "
        println(io, "  ", marker, rpad(string(forcing) * ":", width + 2), " ", origin.label)
    end
end

function _print_overrides(io::IO, data::NamedTuple, indent::String)
    canonical_keys = filter(k -> !_is_special_key(k), keys(data))
    if isempty(canonical_keys)
        return false
    end
    println(io, indent, "overrides (in data):")
    for k in canonical_keys
        println(io, indent, "  ", k)
    end
    return true
end

function _print_problem_init(io::IO, problem)
    init = problem.init
    init === nothing && return
    println(io)
    println(io, "init:")
    for k in keys(init)
        v = init[k]
        summary_str = v === nothing ? "nothing (solver default)" :
                      (v isa AbstractArray ? "$(eltype(v))[$(length(v))]" : repr(v))
        println(io, "  ", rpad(string(k) * ":", 20), " ", summary_str)
    end
end
