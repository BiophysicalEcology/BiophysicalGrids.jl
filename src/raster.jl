# micro_raster.jl
#
# Grid-mode microclimate driver. `MicroRasterProblem` pairs a `MicroMapModel`
# with a rectangular spatial extent, a year range, and an explicit
# `template` grid (a Raster or an RDS source loaded over `area`); `init`
# loads forcings at native resolution, resamples weather and DEM-derived
# terrain (slope, aspect, horizon angles) onto the template, and hands
# off to the shared scaffolding in `micro_common.jl`.
#
# `data::NamedTuple` on the problem provides two layers of override:
#   * collection-level — `weather::RasterStack`, `landcover::RasterStack|Raster`,
#     replacing the corresponding source loader output entirely
#   * single-layer — `dem`, `surface_albedo`, `roughness_height`, or any
#     canonical weather variable name (e.g. `vapour_pressure_deficit`,
#     `mean_temperature`), as a Raster in canonical units; resampled to the
#     run template at `init` time

"""
    MicroRasterProblem(; model, area, dates, template, soil_profile, init=nothing, data=(;))

Concrete grid-microclimate run: pairs a `MicroMapModel` with a spatial
extent, date range, soil column profile, initial conditions, and any
per-run data overrides.

- `model::MicroMapModel`
- `area` — `Extents.Extent` (run every pixel of the template) or any
  GeoInterface-conformant geometry (rasterised into a Bool mask via
  `Rasters.boolmask`; pixels outside the geometry are skipped and left
  as `missing` in the output)
- `dates` — any contiguous date range:
    * `Date(2000, 6, 29)` — single day
    * `Date(2000, 6, 1):Day(1):Date(2000, 6, 30)` — one month
    * `Date(2000, 1, 1):Day(1):Date(2000, 12, 31)` — full year
  Feb 29 is simulated like any other calendar day. Cross-year ranges are supported.
- `template` — spatial grid the run executes on. Required; no fallback.
    * `Type{<:RasterDataSource}` (e.g. `SRTM`) — load that dataset over
      `area` and use it as the grid; weather is resampled onto it.
    * `Raster` — use the supplied raster's `X`/`Y` lookup as the grid;
      weather is resampled onto it. Only `X`/`Y` dims are read.
- `soil_profile::SoilProfile` — per-depth `bulk_density` and
  `mineral_density`. Currently uniform across pixels (a per-pixel
  `soil_profile_source` is a planned extension).
- `init` — initial conditions; default fills `soil_moisture` from
  `model.micro_model.depths` and lets `soil_temperature` fall back to the
  day-mean reference air temperature
- `data::NamedTuple` — overrides for any of:
    * collection sources: `weather` (RasterStack), `landcover` (Stack/Raster)
    * single layers: `dem`, `surface_albedo`, `roughness_height`, or any
      canonical weather-variable name (e.g. `vapour_pressure_deficit`,
      `mean_temperature`, `cloud_cover`). Each as a `Raster` in canonical
      units; resampled to the run template automatically.
    * `shade` — fraction (0–1) of shade cast over the site; default 0.
"""
@kwdef struct MicroRasterProblem{M<:MicroMapModel,A,DT<:Union{Date,AbstractRange{Date}},T,SP<:SoilProfile,IT,D<:NamedTuple,TC<:Timestep}
    model::M
    area::A
    dates::DT
    template::T
    soil_profile::SP
    init::IT = nothing
    data::D = (;)
    # Target output timestep — the step within a day. Default hourly.
    timestep::TC = Hourly()
end

# Positional-`model` convenience: `MicroRasterProblem(model; area, dates, ...)`.
MicroRasterProblem(model::MicroMapModel; kwargs...) =
    MicroRasterProblem(; model, kwargs...)

# ---------------------------------------------------------------------------
# Grid-mode spatial extractors
# ---------------------------------------------------------------------------


# Resample a native-resolution surface property onto the grid template.
# Scalar / Raster cases are handled here; the landcover-weighted path runs
# `_resolve_surface_native` first to compute a native Raster, then passes
# it through this same extractor.
_to_grid_template(value, template) =
    Raster(Fill(value, size(template, X), size(template, Y)),
           dims(template, (X, Y)))
_to_grid_template(r::Raster, template) = Rasters.resample(r; to = template)
@inline function _to_grid_template(r::AbstractArray{<:Quantity}, template)
    u = unit(eltype(r))
    stripped = ustrip.(u, r)
    resampled = Rasters.resample(rebuild(r, stripped); to = template)
    return rebuild(resampled, parent(resampled) .* u)
end

function _resolve_surface_grid(data_override, source, landcover_source, template, area, default_fn)
    native = data_override === nothing ?
        _resolve_surface_native(source, landcover_source, area, default_fn) :
        data_override
    return _to_grid_template(native, template)
end

# Pre-resample every canonical-override raster onto the run template once,
# so the per-pixel hot path is just an `Ti(k)` slice fetch.
function _resample_canonical_overrides(canonical_data::NamedTuple, template)
    return map(canonical_data) do raster
        _to_grid_template(raster, template)
    end
end

# Resolve the user-supplied `MicroRasterProblem.template` into a 2-D Raster.
#   Raster              → use directly (slice off `Ti` if it carries one).
#   Type{<:RDS source}  → load over `area` and use that as the template.
_resolve_template(r::Raster, _area) = hasdim(r, Ti) ? view(r, Ti(1)) : r
# Only the X/Y dims are read downstream (via `Rasters.resample(; to=template)`
# and `dims(template, (X, Y))`), so keep the raster lazy — never materialise.
# Pass `extent` only if the source advertises it in `getraster_keywords`
# (e.g. SRTM mosaics tile selection on it); otherwise load the whole
# dataset lazily and let `crop` trim down to `area`.
function _resolve_template(source::Type{<:RasterDataSources.RasterDataSource}, area::Extent)
    kwargs = :extent in RasterDataSources.getraster_keywords(source) ?
        (; extent = area, lazy = true) : (; lazy = true)
    cropped = crop(Raster(source; kwargs...); to = area, touches = true)
    # Pass to _resolve_template in case it needs Ti dim removed
    return _resolve_template(cropped, nothing)
end

# Degree buffer applied to the weather load area. Needs to be wide enough
# that the cubic-spline resample kernel (4-cell footprint) has source
# cells beyond the requested template extent.
#
# The default (0.5°) covers 1/24° (TerraClimate), 1/120° (CHELSA,
# WorldClim), and 0.25° (ERA5) native grids. Coarser sources override
# this — NCEP's T62 Gaussian grid is ~1.9°, so 4° gives the spline two
# source cells beyond the template extent in every direction.
weather_area_buffer(::Type{<:RasterDataSources.RasterDataSource}) = 0.5
weather_area_buffer(::Type{<:NCEP}) = 4.0

# Resample every layer of a 3-D `(X, Y, Ti)` weather stack onto a 2-D
# `(X, Y)` template — used when the run grid is finer than the native
# weather grid (e.g. SRTM-resolution template over TerraClimate weather).
# `:cubicspline` interpolates continuously across native cell boundaries
# so the upsampled weather varies smoothly with terrain rather than
# stepping in blocks aligned to the coarse weather grid (GDAL's default
# is nearest, which would produce visible block artifacts).
function _resample_weather_to_template(weather::RasterStack, template_2d)
    # GDAL requires (1) a CRS and (2) Sampled{Regular, Intervals(Center())} lookups
    # to compute a valid affine geotransform. Both conditions may be missing for
    # NetCDF-derived data (e.g. CRUCL2 has no embedded CRS; NCEP has Irregular Y).
    # Regularise both template and each source layer before warping.
    warp_crs     = something(crs(template_2d), EPSG(4326))
    eff_template = _regularise_for_resample(template_2d, warp_crs)
    names = keys(weather)
    resampled = map(names) do n
        @info "  resampling layer :$n"
        layer = _regularise_for_resample(getproperty(weather, n), warp_crs)
        read(Rasters.resample(layer; to = eff_template, method = :cubicspline))
    end
    return RasterStack(NamedTuple{names}(resampled))
end

# Rebuild X and Y to Sampled{Regular, Intervals(Center())} and stamp `target_crs`
# so GDAL can always compute a valid affine geotransform.  Applied to both source
# and template layers in `_resample_weather_to_template`.  The rebuild is metadata-
# only (no data copy); `_regular_sampled` reconstructs the same spacing from the
# existing coordinate values.
function _regularise_for_resample(layer::AbstractRaster, target_crs)
    new_x = _regular_sampled(lookup(layer, X))
    new_y = _regular_sampled(lookup(layer, Y))
    return setcrs(Rasters.set(layer, X => new_x, Y => new_y), target_crs)
end

@inline _is_irregular(lk) = span(lk) isa Irregular

# Build a `Sampled{Regular}` lookup whose values land at the mean spacing
# between the existing coordinates. Reversed inputs (descending Y) keep
# their order. `Intervals(Center)` matches GDAL's cell-centre convention.
function _regular_sampled(lk)
    vals = collect(lk)
    n = length(vals)
    if n < 2
        return Sampled(Float64.(vals); sampling = Intervals(Center()),
            order = ForwardOrdered(), span = Regular(1.0))
    end
    step = (Float64(vals[end]) - Float64(vals[1])) / (n - 1)
    rng = range(Float64(vals[1]); step = step, length = n)
    ord = step < 0 ? ReverseOrdered() : ForwardOrdered()
    return Sampled(rng; sampling = Intervals(Center()),
        order = ord, span = Regular(step))
end

# ---------------------------------------------------------------------------
# CommonSolve.init
# ---------------------------------------------------------------------------

"""
    CommonSolve.init(problem::MicroRasterProblem) -> MicroMapCache

Load all per-run data (weather, DEM, landcover, surface properties, canonical
overrides), pre-resolve every grid onto the run template, allocate the
worker-cache pool, and prepare for `solve!`.
"""
function CommonSolve.init(problem::MicroRasterProblem)
    (; model, area, dates, soil_profile, data) = problem
    (; dem_source, weather_source, landcover_source,
       surface_albedo_source, roughness_height_source, soil_moisture_source) = model

    # Inject soil_moisture_source into data before canonical-override resolution,
    # unless the user already supplied data.soil_moisture explicitly.
    # Skipped in solar-only mode — weather is not loaded at all.
    if !model.solar_only && soil_moisture_source !== nothing && !haskey(data, :soil_moisture)
        data = merge(data, (; soil_moisture = _load_prescribed(
            soil_moisture_source, SoilMoisture(), _to_extent(area), _years_from_dates(dates))))
    end

    init_inputs = _resolve_init(problem.init, model.micro_model)
    soil_moisture_available = model.solar_only ? false :
        _has_canonical_input(:soil_moisture, weather_source, data)

    # Normalise the user-supplied dates to a sorted Vector{Date}.
    dates_vec = _normalise_dates(dates)
    years = _years_from_dates(dates)
    calendar = weather_calendar(weather_source)
    native = native_timestep(weather_source)
    target = problem.timestep
    # Ti slicing is in native-timestep units (how the loaded stack is laid out).
    ti_start, ti_end, days_doy = _ti_range_for_dates(calendar, native, years, dates_vec)
    anchor_dates = _step_anchor_dates(calendar, dates_vec)

    # `area` may be a geometry; loaders need an Extent, but the mask uses the geometry itself.
    extent = _to_extent(area)
    @info "init: weather source:   $(weather_source) ($(nameof(typeof(calendar))), $(nameof(typeof(native))) → $(nameof(typeof(target))))"
    @info "init: DEM source:       $(haskey(data, :terrain) ? "skipped (terrain override)" : string(dem_source))"
    @info "init: surface albedo:   $(_source_label(surface_albedo_source))"
    @info "init: roughness height: $(_source_label(roughness_height_source))"

    template_2d = _resolve_template(problem.template, extent)

    # Terrain: use the pre-computed override when supplied (skips DEM download
    # and the expensive horizon-angle sweep). Compute from scratch otherwise.
    terrain = if haskey(data, :terrain)
        @info "init: terrain:          reusing pre-computed terrain from data override"
        data.terrain
    else
        @info "init: terrain:          $(model.compute_terrain ? "slope/aspect/horizon computed from DEM" : "flat (compute_terrain=false)")"
        @info "init: loading DEM..."
        dem_native = _resolve_dem(data, dem_source, extent)
        @info "init: building terrain..."
        if model.compute_terrain
            compute_terrain_grids(dem_native; template = template_2d, n_horizon_angles = N_HORIZON_ANGLES)
        else
            # No terrain computation: crop DEM to the template extent so terrain
            # dims exactly match the run grid (DEM was loaded with a buffer to
            # support slope/aspect at boundaries — not needed here).
            dem_at_template = Rasters.resample(dem_native; to = template_2d, method = :average)
            _flat_terrain_stack(dem_at_template, N_HORIZON_ANGLES)
        end
    end

    # In clear-sky solar-only mode skip weather loading entirely — terrain drives the
    # extent and missing-pixel check. When cloud_correct_solar=true, weather is needed
    # for the cloud cover field. In full mode, load, resample, and slice weather.
    weather = if model.solar_only && !model.cloud_correct_solar
        @info "init: solar-only mode — weather skipped"
        RasterStack()
    else
        buffer_deg = weather_area_buffer(weather_source)
        weather_area = Extents.buffer(extent, (X = buffer_deg, Y = buffer_deg))
        @info "init: loading weather..."
        wf = _resolve_weather(data, weather_source, weather_area, years)
        @info "init: resampling weather to template..."
        wf = _resample_weather_to_template(wf, template_2d)
        @info "init: slicing weather to date range..."
        # Positional integer Ti indexing works for both integer Ti (most sources)
        # and DateTime Ti (ERA5 Zarr). `stack[Ti(...)]` applies the slice to all
        # layers simultaneously — unlike `map(f, stack)` which is pixel-wise.
        wf[Ti(ti_start:ti_end)]
    end

    # Surface-property resampling target: the sliced weather template in full mode
    # (has the right grid after resample); template_2d in clear-sky solar-only mode.
    template = (model.solar_only && !model.cloud_correct_solar) ? template_2d : first(values(weather))
    @info "init: building surface grids..."
    albedo_data = get(data, :surface_albedo, nothing)
    roughness_data = get(data, :roughness_height, nothing)
    albedo_grid = _resolve_surface_grid(albedo_data, surface_albedo_source,
        landcover_source, template, extent, default_landcover_albedo)
    roughness_grid = _resolve_surface_grid(roughness_data, roughness_height_source,
        landcover_source, template, extent, default_landcover_roughness)

    canonical_overrides = _resample_canonical_overrides(_canonical_data(data), template)

    @info "init: building active mask..."
    mask = boolmask(terrain.elevation)
    for l in values(weather)
        mask!(mask; with = hasdim(l, Ti) ? view(l, Ti(1)) : l)
    end
    area isa Extent || mask!(mask; with = area)
    for l in values(canonical_overrides)
        mask!(mask; with = hasdim(l, Ti) ? view(l, Ti(1)) : l)
    end
    @info "init: pixels:          $(count(mask)) active / $(length(mask)) total"

    cloud_constants = _build_cloud_constants()
    solar_pairs = _make_solar_pairs(
        _effective_solar_layers(model), cloud_constants.solar_model.wavelengths)
    @info "init: building cache pool..."
    build_inputs, cache_pool = _build_inputs_and_pool(;
        model, weather_source, weather, terrain, mask,
        albedo_grid, roughness_grid, canonical_overrides,
        init_inputs, soil_moisture_available, years, days = days_doy, cloud_constants,
        soil_profile, target_timestep = target,
    )
    @info "init: done"

    return MicroMapCache(
        problem, weather, terrain, albedo_grid, roughness_grid,
        canonical_overrides, mask, cache_pool,
        (; init_inputs, build_inputs, anchor_dates, solar_pairs),
        cloud_constants,
    )
end

"""
    CommonSolve.solve(problem::MicroRasterProblem) -> RasterStack

Shortcut for `solve!(init(problem))`.
"""
CommonSolve.solve(problem::MicroRasterProblem) = solve!(init(problem))

# ---------------------------------------------------------------------------
# Show methods
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", problem::MicroRasterProblem)
    println(io, "MicroRasterProblem")
    println(io, "  area:     ", problem.area)
    println(io, "  dates:    ", problem.dates)
    println(io, "  template: ", _template_label(problem.template))
    println(io)
    println(io, "forcings:")
    _print_forcing_table(io, _forcing_origins(problem))
    _print_overrides(io, problem.data, "")
    _print_problem_init(io, problem)
end

function Base.show(io::IO, ::MIME"text/plain", cache::MicroMapCache{<:MicroRasterProblem})
    problem = cache.problem
    nx, ny = size(cache.terrain.elevation)
    println(io, "MicroMapCache  (", nx, "×", ny, " grid, ", nx * ny, " pixels)")
    println(io, "  cache_pool:  ", cache.cache_pool.sz_max, " workers")
    println(io)
    println(io, "forcings:")
    _print_forcing_table(io, _forcing_origins(problem))
    _print_overrides(io, problem.data, "")
end


# One-line description of the template for `show`.
_template_label(r::Raster) = "user Raster ($(size(r, X))×$(size(r, Y)))"
_template_label(source::Type{<:RasterDataSources.RasterDataSource}) = string(source)
