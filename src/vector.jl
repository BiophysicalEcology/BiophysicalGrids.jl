# micro_vector.jl
#
# Points-mode microclimate driver. `MicroVectorProblem` pairs a
# `MicroMapModel` with an explicit list of (x, y) coordinates and a year
# range; `init` computes the bounding extent from the points, loads weather
# and DEM at their native resolutions (no resampling), computes terrain at
# native DEM resolution, and per-point extracts every forcing into a 1-D
# `Raster` whose only spatial dim is `Dim{:point}` backed by a
# `MergedLookup` carrying the point coordinates.
#
# Points are taken as any `AbstractVector` of GeoInterface-conformant point
# objects (e.g. `(lon, lat)` tuples, `Point((x, y))`, `(X=…, Y=…)`
# NamedTuples). Coordinate access goes through `GeoInterface.x` / `.y`.

const GI = GeoInterface

# Degree buffer applied to the points' bounding box before passing it to
# weather/landcover loaders. Has to be larger than the coarsest weather
# source's native cell size so `crop(...; touches=true)` returns ≥1 cell
# on Points-sampled grids. 2° accommodates NCEP's ~1.9° Gaussian spacing
# with margin; finer sources just load a few extra cells (per-point Near
# extraction picks the right one).
const _POINTS_LOAD_BUFFER = 2.0

"""
    MicroVectorProblem(; model, points, dates, soil_profile, init=nothing, data=(;))

Concrete points-mode microclimate run: pairs a `MicroMapModel` with a list
of GeoInterface-conformant points, a date range, soil column profile,
initial conditions, and any per-run data overrides.

- `model::MicroMapModel`
- `points::AbstractVector` — any iterable of GeoInterface point-like objects
  (e.g. `Vector{Tuple{Float64,Float64}}` of `(longitude, latitude)`).
  The bounding extent is derived via `GeoInterface.extent(MultiPoint(points))`.
- `dates` — any contiguous date range:
    * `Date(2000, 6, 29)` — single day
    * `Date(2000, 6, 1):Day(1):Date(2000, 6, 30)` — one month
    * `Date(2000, 1, 1):Day(1):Date(2000, 12, 31)` — full year
  Feb 29 is simulated like any other calendar day. Cross-year ranges are supported.
- `soil_profile::SoilProfile` — per-depth `bulk_density` and
  `mineral_density`. Currently uniform across points (a per-point
  `soil_profile_source` is a planned extension).
- `init`, `data` — same semantics as `MicroRasterProblem`. Override Rasters
  (`data.surface_albedo`, `data.roughness_height`, any canonical weather
  variable) are per-point `Near`-extracted at `init` time rather than
  resampled.

The output `RasterStack` is shaped `(Dim{:point}, Ti, [depth | height])`
where the `:point` dim carries a `MergedLookup` of `(x, y)` tuples so
callers can recover coordinates via `lookup`/`val` or by indexing with
`X(At(x)), Y(At(y))` selectors.
"""
@kwdef struct MicroVectorProblem{M<:MicroMapModel,PT<:AbstractVector,DT<:Union{Date,AbstractRange{Date}},SP<:SoilProfile,IT,D<:NamedTuple,TC<:Timestep}
    model::M
    points::PT
    dates::DT
    soil_profile::SP
    init::IT = nothing
    data::D = (;)
    # Target output timestep — the step within a day. Default hourly.
    timestep::TC = Hourly()
end

# Positional-`model` convenience: `MicroVectorProblem(model; points, dates, ...)`.
MicroVectorProblem(model::MicroMapModel; kwargs...) =
    MicroVectorProblem(; model, kwargs...)

# ---------------------------------------------------------------------------
# Points dim + per-point extractors
# ---------------------------------------------------------------------------

# Build the `Dim{:point}` lookup from the user's points. The `MergedLookup`
# stores one `(x, y)` Float64 tuple per point so downstream selector
# indexing (`X(At(x)), Y(At(y))`) can recover the original coordinate.
function _make_points_dim(points)
    coords = [(Float64(GI.x(p)), Float64(GI.y(p))) for p in points]
    return Dim{:point}(MergedLookup(coords, (X(), Y())))
end

# Compute the loading extent from the points via GeoInterface — wrap as a
# MultiPoint and ask GI for its extent (falls back to per-point extrema
# when the wrapper has no stored extent).
function _points_extent(points)
    return GI.extent(GI.Wrappers.MultiPoint(collect(points)); fallback = true)
end

# Per-point extractor: returns a Raster whose only spatial dim is
# `points_dim`. Three input forms:
#   - scalar  → 1-D Fill-Raster (no time)
#   - 2-D `(X, Y)` Raster → 1-D `(Dim{:point},)` Raster of per-point values
#   - 3-D `(X, Y, Ti)` Raster → 2-D `(Dim{:point}, Ti)` Raster
# Unit-tagged inputs flow through untouched — extraction is just indexing,
# so the GDAL unit-stripping that grid-mode resampling needs is not required.

_to_points(value, points_dim) =
    Raster(Fill(value, length(points_dim)), (points_dim,))

function _to_points(r::Raster, points_dim)
    if hasdim(r, Ti)
        return _to_points_with_time(r, points_dim)
    else
        return _to_points_constant(r, points_dim)
    end
end

function _to_points_constant(r::Raster, points_dim)
    coords = lookup(points_dim)
    vals = [r[X(Near(x)), Y(Near(y))] for (x, y) in coords]
    return Raster(vals, (points_dim,))
end

function _to_points_with_time(r::Raster, points_dim)
    coords = lookup(points_dim)
    ti_dim = dims(r, Ti)
    n_pt = length(coords)
    n_ti = length(ti_dim)
    T = eltype(r)
    out = Matrix{T}(undef, n_pt, n_ti)
    @inbounds for (i, (x, y)) in enumerate(coords)
        for t in 1:n_ti
            out[i, t] = r[X(Near(x)), Y(Near(y)), Ti(t)]
        end
    end
    return Raster(out, (points_dim, ti_dim))
end

# Per-layer extract from a RasterStack (weather or terrain). Each layer
# goes through `_to_points` so 2-D layers become 1-D and 3-D layers become
# 2-D, all sharing the same `points_dim` spatial axis. `map(f, ::RasterStack)`
# is pixel-wise, so we iterate layer names explicitly to keep this layer-wise.
function _stack_to_points(stack::RasterStack, points_dim)
    names = keys(stack)
    extracted = map(n -> _to_points(getproperty(stack, n), points_dim), names)
    return RasterStack(NamedTuple{names}(extracted))
end

# Used by weather.jl's points-native fast path (_load_field_at_points).

function _extract_lazy_at_points(lazy_r::Raster, points_dim)
    coords = lookup(points_dim)
    open(lazy_r) do orast
        if ndims(orast) >= 3
            other = first(otherdims(orast, (X, Y)))
            cols = [orast[X(Near(x)), Y(Near(y))] for (x, y) in coords]
            Raster(permutedims(reduce(hcat, map(parent, cols))), (points_dim, other))
        else
            vals = [orast[X(Near(x)), Y(Near(y))] for (x, y) in coords]
            Raster(collect(vals), (points_dim,))
        end
    end
end

# Place a per-point daily total as a single event at each day's first hour
# (midnight); other hours stay zero. Used by BARRA-style sources whose
# rainfall arrives daily while the rest of the stack is hourly.
function _daily_to_hourly_midnight_points(daily::Raster, ref::Raster)
    ti = dims(ref, Ti)
    pdim = dims(daily, 1)
    out = zeros(eltype(daily), length(pdim), length(ti))
    for d in axes(daily, 2)
        out[:, (d - 1) * 24 + 1] .= @view daily[:, d]
    end
    return Raster(out, (pdim, ti))
end

# Broadcast a per-point static value across the Ti axis of `ref`.
function _static_to_ti_points(value::Raster, ref::Raster)
    ti = dims(ref, Ti)
    pdim = dims(value, 1)
    return Raster(repeat(parent(value), 1, length(ti)), (pdim, ti))
end

# Resolve a surface property (albedo / roughness) into a 1-D points-mode
# Raster. Same input forms as grid mode — the native-resolution result is
# computed via `_resolve_surface_native` then `_to_points`-extracted.
function _resolve_surface_points(data_override, source, landcover_source,
                                 points_dim, area, default_fn)
    native = data_override === nothing ?
        _resolve_surface_native(source, landcover_source, area, default_fn) :
        data_override
    return _to_points(native, points_dim)
end

# Pre-extract every canonical-override raster at the point coordinates so
# the per-pixel hot path is just `raster[Dim{:point}(p), Ti(k)]`.
function _extract_canonical_overrides(canonical_data::NamedTuple, points_dim)
    return map(canonical_data) do raster
        _to_points(raster, points_dim)
    end
end

# ---------------------------------------------------------------------------
# CommonSolve.init
# ---------------------------------------------------------------------------

"""
    CommonSolve.init(problem::MicroVectorProblem) -> MicroMapCache

Compute the bounding extent of `problem.points`, load weather and DEM at
native resolution, compute terrain at native DEM resolution, per-point
extract every forcing into `Dim{:point}`-keyed Rasters, allocate the
worker-cache pool, and prepare for `solve!`. No resampling is performed.
"""
function CommonSolve.init(problem::MicroVectorProblem)
    (; model, points, dates, soil_profile, data) = problem
    (; dem_source, weather_source, landcover_source,
       surface_albedo_source, roughness_height_source, soil_moisture_source, init_source) = model

    dates_vec = _normalise_dates(dates)
    years = _years_from_dates(dates)

    # Inject soil_moisture_source into data before canonical-override resolution,
    # unless the user already supplied data.soil_moisture explicitly.
    # Skipped in solar-only mode — weather is not loaded (or is irrelevant).
    if !model.solar_only && soil_moisture_source !== nothing && !haskey(data, :soil_moisture)
        sm_buffer = points_load_buffer(soil_moisture_source)
        sm_area = Extents.buffer(_points_extent(points), (X = sm_buffer, Y = sm_buffer))
        data = merge(data, (; soil_moisture = _load_prescribed(soil_moisture_source, SoilMoisture(), sm_area, years)))
    end

    init_inputs = _resolve_init(problem.init, model.micro_model)
    soil_moisture_available = model.solar_only ? false :
        _has_canonical_input(:soil_moisture, weather_source, data)
    calendar = weather_calendar(weather_source)
    native = native_timestep(weather_source)
    target = problem.timestep
    ti_start, ti_end, days_doy = _ti_range_for_dates(calendar, native, years, dates_vec)
    anchor_dates = _step_anchor_dates(calendar, dates_vec)

    @info "init: weather source:   $(weather_source) ($(nameof(typeof(calendar))), $(nameof(typeof(native))) → $(nameof(typeof(target))))"
    @info "init: DEM source:       $(dem_source)"
    @info "init: surface albedo:   $(_source_label(surface_albedo_source))"
    @info "init: roughness height: $(_source_label(roughness_height_source))"
    @info "init: terrain:          $(model.compute_terrain ? "slope/aspect/horizon computed from DEM" : "flat (compute_terrain=false)")"
    @info "init: points:           $(length(points))"

    area = _points_extent(points)
    points_dim = _make_points_dim(points)

    init_buffer = points_load_buffer(init_source)
    init_soil_native = model.solar_only ? (; soil_temperature = nothing, soil_moisture = nothing) :
        _load_init_soil(init_source,
            Extents.buffer(area, (X = init_buffer, Y = init_buffer)), first(dates_vec))
    init_soil = (;
        soil_temperature = init_soil_native.soil_temperature === nothing ?
            nothing : _to_points(init_soil_native.soil_temperature, points_dim),
        soil_moisture = init_soil_native.soil_moisture === nothing ?
            nothing : _to_points(init_soil_native.soil_moisture, points_dim),
    )

    # Terrain: reuse pre-computed terrain when `data.terrain` is supplied
    # (skips DEM download and horizon-angle sweep). The override must already
    # be per-point (i.e. from a previous `MicroVectorProblem` cache's `terrain`
    # field) — no per-point extraction is performed on it.
    terrain = if haskey(data, :terrain)
        data.terrain
    else
        dem_native = _resolve_dem(data, dem_source, area)
        terrain_native = model.compute_terrain ?
            compute_terrain_grids(dem_native; template = nothing, n_horizon_angles = N_HORIZON_ANGLES) :
            _flat_terrain_stack(dem_native, N_HORIZON_ANGLES)
        _stack_to_points(terrain_native, points_dim)
    end

    # Weather: load at native resolution over the points' bounding extent,
    # buffered so the crop is non-empty even when (a) the points' bbox is
    # degenerate (single point or co-located points) or (b) the weather
    # source uses `Points`-sampling on a coarse grid where `touches=true`
    # would otherwise return zero cells (e.g. NCEP ~2° Gaussian).

    # Skip entirely in clear-sky solar-only mode; load when cloud_correct_solar
    # needs cloud cover. Terrain must be resolved first (above) regardless.
    weather = if model.solar_only && !model.cloud_correct_solar
        RasterStack()
    elseif !haskey(data, :weather) && supports_points_loading(weather_source)
        weather_native = _load_weather_points(weather_source, points_dim, years)
        weather_native[Ti(ti_start:ti_end)]
    else
        weather_buffer = points_load_buffer(weather_source)
        weather_area = Extents.buffer(area, (X = weather_buffer, Y = weather_buffer))
        weather_native = _resolve_weather(data, weather_source, weather_area, years)
        # Slice to the requested date range before per-point extraction —
        # positional Ti indexing works for both integer and DateTime Ti axes.
        weather_native = weather_native[Ti(ti_start:ti_end)]
        _stack_to_points(weather_native, points_dim)
    end

    albedo_data = get(data, :surface_albedo, nothing)
    roughness_data = get(data, :roughness_height, nothing)
    albedo_grid = _resolve_surface_points(albedo_data, surface_albedo_source,
        landcover_source, points_dim, area, default_landcover_albedo)
    roughness_grid = _resolve_surface_points(roughness_data, roughness_height_source,
        landcover_source, points_dim, area, default_landcover_roughness)

    canonical_overrides = _extract_canonical_overrides(_canonical_data(data), points_dim)

    mask = Raster(trues(length(points_dim)), (points_dim,))

    cloud_constants = _build_cloud_constants()
    solar_pairs = _make_solar_pairs(
        _effective_solar_layers(model), cloud_constants.solar_model.wavelengths)
    build_inputs, cache_pool = _build_inputs_and_pool(;
        model, weather_source, weather, terrain, mask,
        albedo_grid, roughness_grid, canonical_overrides,
        init_inputs, soil_moisture_available, years, days = days_doy, cloud_constants,
        soil_profile, target_timestep = target, init_soil,
    )

    return MicroMapCache(
        problem, weather, terrain, albedo_grid, roughness_grid,
        canonical_overrides, mask, cache_pool,
        (; init_inputs, build_inputs, anchor_dates, solar_pairs),
        cloud_constants,
    )
end

"""
    CommonSolve.solve(problem::MicroVectorProblem) -> RasterStack

Shortcut for `solve!(init(problem))`.
"""
CommonSolve.solve(problem::MicroVectorProblem) = solve!(init(problem))

# ---------------------------------------------------------------------------
# Show methods
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", problem::MicroVectorProblem)
    println(io, "MicroVectorProblem")
    println(io, "  points: ", length(problem.points), " point(s)")
    println(io, "  dates:  ", problem.dates)
    println(io)
    println(io, "forcings:")
    _print_forcing_table(io, _forcing_origins(problem))
    _print_overrides(io, problem.data, "")
    _print_problem_init(io, problem)
end

function Base.show(io::IO, ::MIME"text/plain", cache::MicroMapCache{<:MicroVectorProblem})
    problem = cache.problem
    npts = length(problem.points)
    println(io, "MicroMapCache  (", npts, " point", npts == 1 ? "" : "s", ")")
    println(io, "  cache_pool:  ", cache.cache_pool.sz_max, " workers")
    println(io)
    println(io, "forcings:")
    _print_forcing_table(io, _forcing_origins(problem))
    _print_overrides(io, problem.data, "")
end
