abstract type Aspect end
struct Maximum <: Aspect end
struct Minimum <: Aspect end
struct Mean <: Aspect end
struct DiurnalRange <: Aspect end
struct DailyTotal <: Aspect end

# Local stubs for time-of-day qualifiers. Upstream `TimeOfDay` /
# `Microclimate.ClockTime` types are on an unmerged branch; the refactor
# uses them opaquely, so a minimal local declaration matches until those
# land.
abstract type TimeOfDay end
struct ClockTime <: TimeOfDay
    hour::Float64
end
ClockTime(h::Int) = ClockTime(Float64(h))

const Qualifier = Union{Nothing, Aspect, TimeOfDay}

abstract type Sample end

# Each physical quantity a sample can carry: (type name, physical_quantity key,
# working unit). A `nothing` unit is dimensionless. The struct, its zero-arg
# constructor, `physical_quantity`, and `working_unit` are all generated from
# this one table.
const _SAMPLE_QUANTITIES = (
    (:Temperature, :temperature, u"K"),
    (:WindSpeed, :wind_speed, u"m/s"),
    (:EastwardWindSpeed, :eastward_wind, u"m/s"),
    (:NorthwardWindSpeed, :northward_wind, u"m/s"),
    (:RelativeHumidity, :humidity, nothing),
    (:GlobalRadiation, :global_radiation, u"W/m^2"),
    (:LongwaveRadiation, :longwave_radiation, u"W/m^2"),
    (:Rainfall, :rainfall, u"kg/m^2"),
    (:CloudCover, :cloud_cover, nothing),
    (:DewpointTemperature, :dewpoint_temperature, u"K"),
    (:Pressure, :pressure, u"Pa"),
    (:SpecificHumidity, :specific_humidity, nothing),
    (:VapourPressureDeficit, :vapour_pressure_deficit, u"kPa"),
    (:SoilMoisture, :soil_moisture, nothing),
    (:ActualVapourPressure, :actual_vapour_pressure, u"kPa"),
    (:SeaLevelPressure, :sea_level_pressure, u"Pa"),
    (:SoilTemperature, :soil_temperature, u"K"),
    (:ZenithAngle, :zenith_angle, u"°"),
    (:Elevation, :elevation, u"m"),
)

for (Q, pq, unit) in _SAMPLE_QUANTITIES
    @eval begin
        struct $Q{Q<:Qualifier} <: Sample
            at::Q
        end
        $Q() = $Q(nothing)
        @inline physical_quantity(::$Q) = $(QuoteNode(pq))
        @inline working_unit(::$Q) = $unit
    end
end

struct Reference{S<:Sample} <: Sample
    sample::S
end

@inline qualifier(q::Sample) = q.at
@inline qualifier(r::Reference) = qualifier(r.sample)

@inline working_unit(r::Reference) = working_unit(r.sample)
@inline physical_quantity(r::Reference) = physical_quantity(r.sample)

@inline canonical_name(q::Sample) = _canonical_name(q, qualifier(q))
@inline _canonical_name(q::Sample, ::Nothing) = physical_quantity(q)
@inline _canonical_name(q::Sample, ::Maximum) = Symbol(physical_quantity(q), :_max)
@inline _canonical_name(q::Sample, ::Minimum) = Symbol(physical_quantity(q), :_min)
@inline _canonical_name(q::Sample, ::Mean) = Symbol(:mean_, physical_quantity(q))
@inline _canonical_name(q::Sample, ::DiurnalRange) = Symbol(:diurnal_, physical_quantity(q), :_range)
@inline _canonical_name(q::Sample, ::DailyTotal) = Symbol(physical_quantity(q), :_daily)
@inline _canonical_name(q::Sample, t::TimeOfDay) =
    _timed_name(physical_quantity(q), t)
@inline canonical_name(::SoilTemperature{Mean}) = :deep_soil_temperature
@inline canonical_name(r::Reference) = Symbol(:reference_, canonical_name(r.sample))

Base.@assume_effects :foldable function _timed_name(base::Symbol, t::ClockTime)
    h = t.hour
    tag = isinteger(h) ? string(Int(h)) : replace(string(h), "." => "_")
    return Symbol(base, :_, tag)
end

struct Variable{Name, Field, Q, U, T}
    quantity::Q
    unit::U
    transform::T
end
function Variable(quantity::Sample, field::Symbol, unit = 1, transform = identity)
    name = canonical_name(quantity)
    return Variable{name, field, typeof(quantity), typeof(unit), typeof(transform)}(
        quantity, unit, transform)
end

@inline canonical_name(::Variable{Name}) where {Name} = Name
@inline native_field(::Variable{<:Any, Field}) where {Field} = Field
@inline quantity(v::Variable) = v.quantity

abstract type Calendar end
struct Monthly <: Calendar end
struct Daily <: Calendar end

abstract type Timestep end
struct MinMax <: Timestep end
struct SubDaily{N} <: Timestep end
const Hourly = SubDaily{24}
const SixHourly = SubDaily{4}

@inline samples_per_day(::SubDaily{N}) where {N} = N
@inline samples_per_day(::MinMax) = 1

@inline weather_calendar(::Type) = Monthly()

@inline native_timestep(::Type) = MinMax()

@inline days_per_year(::Monthly) = 12
@inline days_per_year(::Daily) = 365

@inline steps_per_year(cal::Calendar, cad::Timestep) =
    days_per_year(cal) * samples_per_day(cad)

# ---------------------------------------------------------------------------
# Loader traits
# ---------------------------------------------------------------------------
abstract type Loader end

struct YearlyTimeSeries <: Loader end
struct MonthlyTimeSeries <: Loader end
struct MonthlyClimatology <: Loader end
struct MonthlyClimatologyPeriod <: Loader end
struct MultiBandClimatologyPeriod <: Loader end
struct SingleFileBands <: Loader end
struct DailyFiles <: Loader end
struct ContiguousTimeSeries <: Loader end
function loader end

abstract type LongitudeConvention end
struct Longitude180 <: LongitudeConvention end
struct Longitude360 <: LongitudeConvention end

@inline longitude_convention(::Type) = Longitude180()

# Degree buffer for points-mode bounding-box loads (see `_POINTS_LOAD_BUFFER`
# in vector.jl) -- default suits coarse grids (e.g. NCEP ~1.9°). Finer
# sources should override: at ERA5-Land's ~9km spacing, the default buffer
# pulls in a ~40x40 cell block per point instead of a handful.
@inline points_load_buffer(::Type) = _POINTS_LOAD_BUFFER
@inline points_load_buffer(::Nothing) = _POINTS_LOAD_BUFFER  # init_source is optional

@inline _native_lon_crop(::Longitude180, area::Extent) = (area, identity)
function _native_lon_crop(::Longitude360, area::Extent)
    (area.X[1] >= 0 && area.X[2] >= 0) && return (area, identity)
    load_area = Extent(X = (mod(area.X[1], 360.0), mod(area.X[2], 360.0)), Y = area.Y)
    return (load_area, spatial_dims -> _shift_x_lookup(spatial_dims, -360.0))
end

function _shift_x_lookup(spatial_dims, by)
    x_lk = lookup(spatial_dims, X)
    new_x = X(Sampled(collect(x_lk) .+ by;
        order = order(x_lk), span = span(x_lk), sampling = sampling(x_lk)))
    return (new_x, dims(spatial_dims, Y))
end

# ---------------------------------------------------------------------------
# Per-source extension points (with defaults)
# ---------------------------------------------------------------------------

@inline layers(source) = map(canonical_name, variables(source))
# Variables usable for init-only seeding (e.g. ERA5's soil layers, not part
# of its regular per-timestep `variables()`). Defaults to `variables(source)`.
@inline init_variables(source) = variables(source)
@inline fallback_source(::Type) = nothing
# A source's fallback layers are exactly the baseline's layers it does not
# provide itself, so they follow from `variables` and `fallback_source`.
@inline fallback_layers(source) = _fallback_layers(source, fallback_source(source))
@inline _fallback_layers(_source, ::Nothing) = ()
@inline _fallback_layers(source, baseline) = _names_absent(layers(source), layers(baseline))
@inline _names_absent(_own, ::Tuple{}) = ()
@inline _names_absent(own, candidates::Tuple) =
    first(candidates) in own ?
        _names_absent(own, Base.tail(candidates)) :
        (first(candidates), _names_absent(own, Base.tail(candidates))...)
@inline _extra_getraster_kwargs(::Type) = (;)
@inline _extra_getpoint_kwargs(::Type) = (;)

# variables(source) alone misses fallback-contributed quantities (e.g. SILO's
# CRUCL2 wind), so they'd never get copied into buffers.
@inline _full_variables(source) = _merge_fallback_variables(source, fallback_source(source))
@inline _merge_fallback_variables(source, ::Nothing) = variables(source)
@inline _merge_fallback_variables(source, baseline) =
    (variables(source)..., map(name -> _variable_for(variables(baseline), name), fallback_layers(source))...)

@inline supports_points_loading(::Type) = true

# A source's grid elevation is just its declared `Elevation()` quantity —
# static, no Ti axis. Sources that don't declare it get `nothing` and fall
# back to the site DEM elevation.
@inline weather_grid_elevation(source::Type, weather, I) =
    _grid_elevation(_maybe_variable(variables(source), :elevation), weather, I)
@inline _grid_elevation(::Nothing, _weather, _I) = nothing
@inline _grid_elevation(var::Variable, weather, I) =
    u"m"(var.transform(weather[:elevation][I...]) * var.unit)

function _aggregate_ti_to_daily(layer::AbstractRaster, factor::Int)
    return Rasters.aggregate(mean, layer, (Ti(factor),))
end

function _load_weather(source::Type, area::Extent, years)
    primary_stack = _load_canonical(source, layers(source), area, years)
    baseline = fallback_source(source)
    stack = baseline === nothing ?
        RasterStack(primary_stack) :
        RasterStack(merge(primary_stack,
            _match_fallback_resolution(source, baseline,
                _load_canonical(baseline, fallback_layers(source), area, years),
                primary_stack, years)))
    # Source files declare a `missingval`, so loaded eltypes are
    # `Union{Missing, T}`. The per-pixel reader does `value * unit`,
    # which throws `convert(Missing, Quantity)` on any masked cell.
    # `replace_missing(..., NaN)` strips the Missing union and lets
    # NaN propagate visibly if any cell is genuinely masked.
    stack = Rasters.replace_missing(stack, NaN)
    # Ti-count check runs on the first Ti-varying layer only — static ones
    # (Elevation() declarations) have no Ti axis and are excluded.
    _check_ti_length(source, stack, years)
    return stack
end

# A trailing, still-in-progress calendar year (e.g. SILO's near-real-time
# data) legitimately loads fewer days than a full year; only that shortfall
# is allowed through.
function _check_ti_length(source::Type, stack::RasterStack, years)
    cal = weather_calendar(source)
    native = native_timestep(source)
    expected = length(_days_of_year(cal, years)) * samples_per_day(native)
    ref = _first_ti_layer(stack)
    n = length(dims(ref, Ti))
    n == expected && return nothing
    prior_expected = length(_days_of_year(cal, years[1:end-1])) * samples_per_day(native)
    is_trailing_partial_year = cal isa Daily &&
        last(years) == Dates.year(Dates.today()) &&
        prior_expected <= n < expected
    is_trailing_partial_year && return nothing
    error("$(nameof(source)): loaded Ti=$n, but the source declares $expected " *
        "($(nameof(typeof(cal))) × " *
        "$(nameof(typeof(native))), $(length(years)) years)")
end

function _load_weather_points(source::Type, points_dim, years)
    primary_stack = _load_canonical_points(source, layers(source), points_dim, years)
    baseline = fallback_source(source)
    merged = baseline === nothing ?
        primary_stack :
        merge(primary_stack,
            _match_fallback_resolution_points(source, baseline,
                _load_canonical_points(baseline, fallback_layers(source), points_dim, years),
                primary_stack, years))
    stack = RasterStack(_truncate_to_common_ti(merged))
    stack = Rasters.replace_missing(stack, NaN)
    _check_ti_length(source, stack, years)
    return stack
end

@inline _first_ti_layer(stack::RasterStack) =
    first(l for l in values(stack) if hasdim(l, Ti))

# Near-real-time layers (e.g. SILO's current year) can differ in day count
# across variables; trim all Ti layers to the shortest so they merge cleanly.
function _truncate_to_common_ti(layers::NamedTuple)
    ti_lengths = [length(dims(l, Ti)) for l in values(layers) if hasdim(l, Ti)]
    isempty(ti_lengths) && return layers
    n = minimum(ti_lengths)
    all(==(n), ti_lengths) && return layers
    return map(layers) do l
        hasdim(l, Ti) && length(dims(l, Ti)) > n ? view(l, Ti(1:n)) : l
    end
end

# ---------------------------------------------------------------------------
# Per-loader file-reading
# ---------------------------------------------------------------------------

function _load_layers(::YearlyTimeSeries, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    layers = map(fields) do name
        per_year = map(years) do yr
            @info "  loading $source $name $yr..."
            path = getraster(source, name; date = Date(yr), extras...)
            read(crop(Raster(path; name, lazy = true); to = area, touches = true))
        end
        cat(per_year...; dims = Ti)
    end
    return NamedTuple{fields}(layers)
end
function _load_layers(::MonthlyTimeSeries, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    dates = _monthly_date_sequence(years)
    layers = map(fields) do name
        per_month = map(dates) do d
            @info "  loading $source $name $(year(d))-$(month(d))..."
            path = getraster(source, name; date = d, extras...)
            read(crop(Raster(path; name, lazy = true); to = area, touches = true))
        end
        cat(per_month...; dims = Ti)
    end
    return NamedTuple{fields}(layers)
end
function _load_layers(::MonthlyClimatology, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    layers = map(fields) do name
        @info "  loading $source $name (monthly climatology)..."
        _build_monthly_climatology(area, nyears, name) do month
            getraster(source, name; month)
        end
    end
    return NamedTuple{fields}(layers)
end
function _load_layers(::MonthlyClimatologyPeriod, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    period_date = Date(first(years))
    layers = map(fields) do name
        @info "  loading $source $name (monthly climatology, $period_date)..."
        _build_monthly_climatology(area, nyears, name) do month
            getraster(source, name; date = period_date, month)
        end
    end
    return NamedTuple{fields}(layers)
end
function _load_layers(::MultiBandClimatologyPeriod, source, fields::Tuple,
                      area::Extent, years)
    nyears = length(years)
    period_date = Date(first(years))
    layers = map(fields) do name
        @info "  loading $source $name ($period_date)..."
        path = getraster(source, name; date = period_date)
        full = read(crop(Raster(path; name, lazy = true); to = area, touches = true))
        # File holds 12 months along the 3rd dim (Band or Ti, depending on
        # the format). Pull out the data array, tile across years, and re-
        # wrap with a proper `Ti` axis.
        data = parent(full)
        size(data, 3) == 12 || error(
            "$source: expected a 12-month multi-band file for $name, got $(size(data, 3)) bands")
        tiled = nyears == 1 ? data : repeat(data; outer = (1, 1, nyears))
        spatial_dims = dims(full)[1:2]
        Raster(tiled, (spatial_dims..., Ti(1:(12 * nyears))); crs = crs(full))
    end
    return NamedTuple{fields}(layers)
end
function _load_layers(::SingleFileBands, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    load_area, restore_dims = _native_lon_crop(longitude_convention(source), area)
    @info "  loading $source $fields..."
    path = getraster(source; _extra_getraster_kwargs(source)...)
    raw = read(crop(RasterStack(path; name = fields, lazy = true);
                    to = load_area, touches = true))
    # Infer native step count from the first 3-D layer (some layers, e.g. CRUCL2 :elv, are static 2-D).
    nsteps = 1
    for name in fields
        ndims(raw[name]) >= 3 && (nsteps = size(parent(raw[name]), 3); break)
    end
    ntotal = nsteps * nyears
    layers = map(fields) do name
        lyr  = raw[name]
        data = parent(lyr)
        tiled = if ndims(data) == 2
            repeat(reshape(data, size(data)..., 1); outer = (1, 1, ntotal))
        else
            nyears == 1 ? data : repeat(data; outer = (1, 1, nyears))
        end
        Raster(tiled, (restore_dims(dims(lyr)[1:2])..., Ti(1:ntotal)); crs = crs(lyr))
    end
    return NamedTuple{fields}(layers)
end
function _load_layers(::DailyFiles, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    dates = _daily_date_sequence(years)
    nsteps = length(dates)
    layers = map(fields) do name
        @info "  loading $source $name ($nsteps daily files, $(first(dates)) to $(last(dates)))..."
        per_day = map(dates) do d
            path = getraster(source, name; date = d, extras...)
            read(crop(Raster(path; name, lazy = true); to = area, touches = true))
        end
        first_day = first(per_day)
        spatial_dims = dims(first_day)
        stacked = cat(map(parent, per_day)...; dims = 3)
        Raster(stacked, (spatial_dims..., Ti(1:nsteps)); crs = crs(first_day))
    end
    return NamedTuple{fields}(layers)
end
# ARCO-ERA5's Zarr dims come back as bare integer NoLookups (no real
# coordinate values), so coordinates are read directly from the store —
# see ext/MicroclimateMapperZarrExt. Dispatches on the cloud source object,
# not its bare `.url`, so an authenticated source gets its own method.
function _contiguous_series_coords end
_contiguous_series_coords(cloud_source) = error(
    "loading a ContiguousTimeSeries source (e.g. ERA5) requires `using ZarrDatasets`",
)

# `RasterStack(url; source=Zarrsource())`'s field discovery misses most of
# this store's arrays, so named variables are opened directly instead —
# see ext/MicroclimateMapperZarrExt.
function _contiguous_series_open end
_contiguous_series_open(cloud_source, long_name) = error(
    "loading a ContiguousTimeSeries source (e.g. ERA5) requires `using ZarrDatasets`",
)

# lon/lat/time index ranges for one store's coords -- shared by the
# single-store and multi-group loaders. `findall`, not a direction-assuming
# findfirst/findlast pair, since stores disagree on ascending vs descending.
function _contiguous_series_indices(source::Type, coords, area::Extent, time_start::DateTime, time_end::DateTime)
    hstart = Dates.value(time_start - coords.epoch) ÷ coords.step_ms
    hend = Dates.value(time_end - coords.epoch) ÷ coords.step_ms
    ti = searchsortedfirst(coords.hours, hstart):searchsortedlast(coords.hours, hend)

    load_area, _ = _native_lon_crop(longitude_convention(source), area)
    lon_idxs = findall(x -> load_area.X[1] <= x <= load_area.X[2], coords.lon)
    lat_idxs = findall(y -> area.Y[1] <= y <= area.Y[2], coords.lat)
    xi, yi = first(lon_idxs):last(lon_idxs), first(lat_idxs):last(lat_idxs)
    return (; xi, yi, ti, xs = coords.lon[xi], ys = coords.lat[yi])
end

function _load_contiguous_series(source, fields::Tuple, area::Extent, time_start::DateTime, time_end::DateTime)
    cloud_source = getraster(source)
    coords = _contiguous_series_coords(cloud_source)
    (; xi, yi, ti, xs, ys) = _contiguous_series_indices(source, coords, area, time_start, time_end)

    layers = map(fields) do name
        @info "  loading $source $name..."
        long_name = layername(source, name)
        raw = _contiguous_series_open(cloud_source, long_name)
        data = raw[xi, yi, ti]
        Raster(data, (X(xs), Y(ys), Ti(1:length(ti))); crs = EPSG(4326), name)
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::ContiguousTimeSeries, source, fields::Tuple, area::Extent, years)
    time_start = DateTime(first(years), 1, 1, 0)
    time_end = DateTime(last(years), 12, 31, 23)
    _load_contiguous_series(source, fields, area, time_start, time_end)
end

# A source whose fields span several separate base stores (e.g. ERA5ECMWFLand's
# per-topic-group Zarr stores) rather than one shared store. `native_group`
# resolves each field to its group; one `getraster` call per distinct group.
struct MultiGroupContiguousTimeSeries <: Loader end
function native_group end

function _load_multi_group_series(source, fields::Tuple, area::Extent, time_start::DateTime, time_end::DateTime)
    groups = unique(map(f -> native_group(source, f), fields))
    group_sources = Dict(g => getraster(source, g) for g in groups)
    layers = map(fields) do name
        @info "  loading $source $name..."
        cloud_source = group_sources[native_group(source, name)]
        coords = _contiguous_series_coords(cloud_source)
        (; xi, yi, ti, xs, ys) = _contiguous_series_indices(source, coords, area, time_start, time_end)
        long_name = layername(source, name)
        raw = _contiguous_series_open(cloud_source, long_name)
        data = raw[xi, yi, ti]
        Raster(data, (X(xs), Y(ys), Ti(1:length(ti))); crs = EPSG(4326), name)
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::MultiGroupContiguousTimeSeries, source, fields::Tuple, area::Extent, years)
    time_start = DateTime(first(years), 1, 1, 0)
    time_end = DateTime(last(years), 12, 31, 23)
    _load_multi_group_series(source, fields, area, time_start, time_end)
end

"""
    prefetch_weather!(source, points, dates; fields = layers(source), batch = Week(1))

Download and cache `source` weather data for `points`' bounding area and
`dates`, one `batch`-sized NetCDF file at a time under
`RasterDataSources.rasterpath(source)/prefetch`. Skips batches already
cached, so a re-run resumes where it left off.
"""
function prefetch_weather!(source::Type, points, dates; fields = layers(source), batch = Week(1))
    buffer = points_load_buffer(source)
    area = Extents.buffer(_points_extent(points), (X = buffer, Y = buffer))
    date_start, date_end = extrema(dates)
    cache_dir = joinpath(RasterDataSources.rasterpath(source), "prefetch")
    mkpath(cache_dir)

    for b0 in date_start:batch:date_end
        b1 = min(b0 + batch - Day(1), date_end)
        path = _prefetch_path(cache_dir, fields, area, b0, b1)
        if isfile(path)
            @info "prefetch: $b0 to $b1 already cached"
            continue
        end
        @info "prefetch: downloading $b0 to $b1..."
        stack = RasterStack(_load_contiguous_series(source, fields, area,
            DateTime(b0), DateTime(b1) + Day(1) - Second(1)))
        write(path, stack)
    end
    return cache_dir
end

_prefetch_path(dir, fields, area, b0, b1) = joinpath(dir,
    "$(join(fields, '-'))_" *
    "$(round(area.X[1]; digits = 2))_$(round(area.X[2]; digits = 2))_" *
    "$(round(area.Y[1]; digits = 2))_$(round(area.Y[2]; digits = 2))_" *
    "$(Dates.format(b0, "yyyymmdd"))_$(Dates.format(b1, "yyyymmdd")).nc")

function _load_field_at_points(source, name, points_dim; kw...)
    lazy_r = Raster(source, name; lazy = true, kw...)
    return _extract_lazy_at_points(lazy_r, points_dim)
end

function _load_layers_at_points(loader::Loader, source, fields::Tuple, points_dim, years)
    area_layers = _load_layers(loader, source, fields, _points_bbox(source, points_dim), years)
    return NamedTuple{fields}(map(l -> _to_points(l, points_dim), values(area_layers)))
end

function _load_layers_at_points(::YearlyTimeSeries, source, fields::Tuple, points_dim, years)
    extras = _extra_getraster_kwargs(source)
    layers = map(fields) do name
        per_year = map(years) do yr
            @info "  loading $source $name $yr (points)..."
            _load_field_at_points(source, name, points_dim; date = Date(yr), extras...)
        end
        cat(per_year...; dims = Ti)
    end
    return NamedTuple{fields}(layers)
end
function _load_layers_at_points(::MonthlyTimeSeries, source, fields::Tuple, points_dim, years)
    extras = _extra_getraster_kwargs(source)
    dates = _monthly_date_sequence(years)
    layers = map(fields) do name
        per_month = map(dates) do d
            @info "  loading $source $name $(year(d))-$(month(d)) (points)..."
            _load_field_at_points(source, name, points_dim; date = d, extras...)
        end
        cat(per_month...; dims = Ti)
    end
    return NamedTuple{fields}(layers)
end
function _load_layers_at_points(::DailyFiles, source, fields::Tuple, points_dim, years)
    extras = _extra_getraster_kwargs(source)
    dates = _daily_date_sequence(years)
    layers = map(fields) do name
        @info "  loading $source $name ($(length(dates)) daily files, $(first(dates)) to $(last(dates)), points)..."
        per_day = map(dates) do d
            _load_field_at_points(source, name, points_dim; date = d, extras...)
        end
        cat(per_day...; dims = Ti)
    end
    return NamedTuple{fields}(layers)
end

function _daily_date_sequence(years)
    dates = Date[]
    for y in years
        for d in Date(y, 1, 1):Day(1):Date(y, 12, 31)
            push!(dates, d)
        end
    end
    return dates
end

function _monthly_date_sequence(years)
    dates = Date[]
    for y in years
        for m in 1:12
            push!(dates, Date(y, m, 1))
        end
    end
    return dates
end

# Aligns a fallback source's layers onto `source`'s own Ti axis when the two
# calendars differ (e.g. CRUCL2's 12-month climatology backfilling SILO's
# daily wind). Same calendar is a no-op; unhandled combinations raise clearly.
_match_fallback_resolution(source::Type, baseline::Type, layers::NamedTuple, template::NamedTuple, years) =
    _match_fallback_resolution(weather_calendar(source), weather_calendar(baseline), layers, template, years)
_match_fallback_resolution(::T, ::T, layers::NamedTuple, _template, _years) where {T <: Calendar} = layers
function _match_fallback_resolution(::Daily, ::Monthly, layers::NamedTuple, template::NamedTuple, years)
    ref = first(values(template))
    return NamedTuple{keys(layers)}(map(l -> _monthly_to_daily(l, ref, years), values(layers)))
end
_match_fallback_resolution(target::Calendar, base::Calendar, ::NamedTuple, _template, _years) =
    error("No calendar-matching path from $(typeof(base)) to $(typeof(target)) for a weather fallback")

# Points-mode counterpart of `_match_fallback_resolution`.
_match_fallback_resolution_points(source::Type, baseline::Type, layers::NamedTuple, template::NamedTuple, years) =
    _match_fallback_resolution_points(weather_calendar(source), weather_calendar(baseline), layers, template, years)
_match_fallback_resolution_points(::T, ::T, layers::NamedTuple, _template, _years) where {T <: Calendar} = layers
function _match_fallback_resolution_points(::Daily, ::Monthly, layers::NamedTuple, template::NamedTuple, years)
    ref = first(values(template))
    return NamedTuple{keys(layers)}(map(l -> _monthly_to_daily_points(l, ref, years), values(layers)))
end
_match_fallback_resolution_points(target::Calendar, base::Calendar, ::NamedTuple, _template, _years) =
    error("No calendar-matching path from $(typeof(base)) to $(typeof(target)) for a weather fallback (points mode)")

# Resample a monthly climatology onto `ref`'s grid, then repeat each of its
# 12*nyears slices across every real day of that month (leap-aware), reusing
# `ref`'s own Ti axis. Setup-time only. `Rasters.resample` rejects a
# degenerate 1x1 target grid (e.g. a `PointQuery`-sourced `ref`), so that
# case selects the nearest source cell directly instead.
function _monthly_to_daily(layer::AbstractRaster, ref::AbstractRaster, years)
    x_lk, y_lk = dims(ref, X), dims(ref, Y)
    resampled = if length(x_lk) == 1 && length(y_lk) == 1
        cell = layer[X(Near(only(x_lk))), Y(Near(only(y_lk)))]
        Raster(reshape(parent(cell), 1, 1, length(cell)), (x_lk, y_lk, dims(layer, Ti)); crs = crs(layer))
    else
        Rasters.resample(layer; to = (x_lk, y_lk))
    end
    data = parent(resampled)
    ti = dims(ref, Ti)
    out = similar(data, size(data, 1), size(data, 2), length(ti))
    years_v = collect(years)
    nout = length(ti)
    d = 0
    for (yi, y) in enumerate(years_v), m in 1:12
        d >= nout && break
        @views month_slice = data[:, :, (yi - 1) * 12 + m]
        for _ in 1:Dates.daysinmonth(y, m)
            d >= nout && break
            d += 1
            out[:, :, d] .= month_slice
        end
    end
    return Raster(out, (dims(resampled)[1:2]..., ti); crs = crs(resampled))
end

# Points-mode `_monthly_to_daily` -- no resample needed, points already sit
# at their own native-grid location.
function _monthly_to_daily_points(layer::AbstractRaster, ref::AbstractRaster, years)
    data = parent(layer)
    pdim = dims(layer, 1)
    ti = dims(ref, Ti)
    out = similar(data, size(data, 1), length(ti))
    years_v = collect(years)
    nout = length(ti)
    d = 0
    for (yi, y) in enumerate(years_v), m in 1:12
        d >= nout && break
        @views month_col = data[:, (yi - 1) * 12 + m]
        for _ in 1:Dates.daysinmonth(y, m)
            d >= nout && break
            d += 1
            out[:, d] .= month_col
        end
    end
    return Raster(out, (pdim, ti))
end

# Resample a daily layer onto `ref`'s grid, then place each day's value at
# its first hour (midnight), zero elsewhere — a single daily event on `ref`'s
# hourly Ti axis. Setup-time only.
function _daily_to_hourly_midnight(layer::AbstractRaster, ref::AbstractRaster)
    resampled = Rasters.resample(layer; to = dims(ref, (X, Y)))
    data = parent(resampled)
    ti = dims(ref, Ti)
    out = zeros(eltype(data), size(data, 1), size(data, 2), length(ti))
    for d in 1:size(data, 3)
        out[:, :, (d - 1) * 24 + 1] .= @view data[:, :, d]
    end
    return Raster(out, (dims(resampled)[1:2]..., ti); crs = crs(resampled))
end

# Resample a static layer onto `ref`'s grid and tile across `ref`'s Ti axis
# (CRUCL2 `:elv`, GRIDMET `:elev` convention). Setup-time only.
function _static_to_ti(layer::AbstractRaster, ref::AbstractRaster)
    resampled = Rasters.resample(layer; to = dims(ref, (X, Y)))
    data = parent(resampled)
    ti = dims(ref, Ti)
    tiled = repeat(reshape(data, size(data)..., 1); outer = (1, 1, length(ti)))
    return Raster(tiled, (dims(resampled)[1:2]..., ti); crs = crs(resampled))
end


function _build_monthly_climatology(path_for, area::Extent, nyears::Int, name)
    monthly = map(1:12) do m
        path = path_for(m)
        read(crop(Raster(path; name, lazy = true); to = area, touches = true))
    end
    first_monthly = first(monthly)
    spatial_dims = dims(first_monthly)
    climatology = cat(map(parent, monthly)...; dims = 3)
    tiled = nyears == 1 ? climatology :
                                  repeat(climatology; outer = (1, 1, nyears))
    return Raster(tiled, (spatial_dims..., Ti(1:(12 * nyears))); crs = crs(first_monthly))
end

@inline _unique_native_fields(vars::Tuple) = _unique_native_fields((), vars)
@inline _unique_native_fields(acc::Tuple, ::Tuple{}) = acc
@inline function _unique_native_fields(acc::Tuple, vars::Tuple)
    f = native_field(first(vars))
    new_acc = _contains(acc, f) ? acc : (acc..., f)
    return _unique_native_fields(new_acc, Base.tail(vars))
end
@inline _contains(::Tuple{}, _) = false
@inline _contains(t::Tuple, f) = first(t) === f || _contains(Base.tail(t), f)

function _load_canonical(source, names::Tuple, area::Extent, years)
    vars = _select_variables(variables(source), names)
    # `Elevation()` declarations are static — one file, no Ti axis. Split
    # them off and load via the one-shot path; the rest go through the
    # source's declared Ti-varying loader.
    ti_vars, static_vars = _split_static(vars)
    ti_stack     = _load_layers(loader(source), source,
                                _unique_native_fields(ti_vars), area, years)
    static_stack = _load_static(source, _unique_native_fields(static_vars), area)
    return _canonical_keyed(vars, merge(ti_stack, static_stack))
end

function _load_canonical_points(source, names::Tuple, points_dim, years)
    vars = _select_variables(variables(source), names)
    ti_vars, static_vars = _split_static(vars)
    ti_stack = _load_layers_at_points(loader(source), source,
                                      _unique_native_fields(ti_vars), points_dim, years)
    static_fields = _unique_native_fields(static_vars)
    static_stack = if isempty(static_fields)
        NamedTuple()
    else
        area = _points_bbox(source, points_dim)
        native = _load_static(source, static_fields, area)
        NamedTuple{static_fields}(map(l -> _to_points(l, points_dim), values(native)))
    end
    return _canonical_keyed(vars, merge(ti_stack, static_stack))
end

function _points_bbox(source, points_dim)
    coords = lookup(points_dim)
    lons = first.(coords); lats = last.(coords)
    buffer = points_load_buffer(source)
    return Extents.buffer(Extent(X = extrema(lons), Y = extrema(lats)), (X = buffer, Y = buffer))
end

# Compile-time partition of a variables tuple into (Ti-varying, static)
# using method dispatch on the Variable's quantity type parameter.
@inline _split_static(vars::Tuple) = _split_static((), (), vars)
@inline _split_static(ti::Tuple, st::Tuple, ::Tuple{}) = (ti, st)
@inline _split_static(ti::Tuple, st::Tuple, vars::Tuple) =
    _split_static_step(ti, st, first(vars), Base.tail(vars))
@inline _split_static_step(ti::Tuple, st::Tuple,
        v::Variable{<:Any, <:Any, <:Elevation}, rest::Tuple) =
    _split_static(ti, (st..., v), rest)
@inline _split_static_step(ti::Tuple, st::Tuple, v::Variable, rest::Tuple) =
    _split_static((ti..., v), st, rest)

# One-shot read for static (no-Ti) fields. Dispatched on the source's
# loader so per-file-organization access patterns are respected: the
# default is `getraster(source, name)` (BARRA-style, one file per field);
# SingleFileBands reads bands from the source's single bundled file.
@inline _load_static(source, fields::Tuple, area::Extent) =
    _load_static(loader(source), source, fields, area)

_load_static(_loader, _source, ::Tuple{}, _area) = NamedTuple()
# Disambiguates against each loader-specific method's ::Tuple{} match below.
_load_static(::Loader, _source, ::Tuple{}, _area::Extent) = NamedTuple()
_load_static(::SingleFileBands, _source, ::Tuple{}, _area::Extent) = NamedTuple()
function _load_static(::Loader, source, fields::Tuple, area::Extent)
    layers = map(fields) do name
        @info "  loading $source $name (static)..."
        path = getraster(source, name)
        read(crop(Raster(path; name, lazy = true); to = area, touches = true))
    end
    return NamedTuple{fields}(layers)
end
function _load_static(::SingleFileBands, source, fields::Tuple, area::Extent)
    load_area, restore_dims = _native_lon_crop(longitude_convention(source), area)
    path = getraster(source; _extra_getraster_kwargs(source)...)
    raw = read(crop(RasterStack(path; name = fields, lazy = true);
                    to = load_area, touches = true))
    layers = map(fields) do name
        lyr = raw[name]
        Raster(parent(lyr), restore_dims(dims(lyr)[1:2]); crs = crs(lyr))
    end
    return NamedTuple{fields}(layers)
end

@inline _select_variables(vars::Tuple, names::Tuple) =
    map(name -> _variable_for(vars, name), names)
@inline _variable_for(::Tuple{}, name) =
    error("no weather variable declares the canonical name :$name")
@inline _variable_for(vars::Tuple, name) =
    canonical_name(first(vars)) === name ? first(vars) :
    _variable_for(Base.tail(vars), name)

@inline _maybe_variable(::Tuple{}, _name) = nothing
@inline _maybe_variable(vars::Tuple, name) =
    canonical_name(first(vars)) === name ? first(vars) :
    _maybe_variable(Base.tail(vars), name)

@inline _canonical_keyed(vars::Tuple, native_stack) =
    NamedTuple{map(canonical_name, vars)}(
        map(v -> getproperty(native_stack, native_field(v)), vars))

function _load_prescribed(source, sample::Sample, area::Extent, years)
    name = canonical_name(sample)
    raw = Rasters.replace_missing(_load_canonical(source, (name,), area, years)[name], NaN)
    var = _variable_for(variables(source), name)
    converted = ustrip.(canonical_unit(name), var.transform.(parent(raw)) .* var.unit)
    return Raster(converted, dims(raw); crs = crs(raw))
end

# Cheap single-snapshot (X, Y) grid near `start_date`, for seeding initial
# conditions from a source that isn't the run's `weather_source` (e.g. ERA5
# seeding an NCEP/AWAP/SILO run). Not a time-varying forcing: a
# ContiguousTimeSeries/MultiGroupContiguousTimeSeries sources are read for a
# single hour; other loaders read a single year.
function _load_init_snapshot(source, sample::Sample, area::Extent, start_date::Date)
    name = canonical_name(sample)
    var = _variable_for(init_variables(source), name)
    field = native_field(var)
    raw2d = if loader(source) isa ContiguousTimeSeries
        t0 = DateTime(start_date)
        first(values(_load_contiguous_series(source, (field,), area, t0, t0 + Hour(1))))[Ti(1)]
    elseif loader(source) isa MultiGroupContiguousTimeSeries
        t0 = DateTime(start_date)
        first(values(_load_multi_group_series(source, (field,), area, t0, t0 + Hour(1))))[Ti(1)]
    else
        yr = year(start_date)
        first(values(_load_layers(loader(source), source, (field,), area, yr:yr)))[Ti(1)]
    end
    raw2d = Rasters.replace_missing(raw2d, NaN)
    # Kept Unitful, unlike _load_prescribed's ustrip -- MicroInputs wants Unitful K.
    converted = var.transform.(parent(raw2d)) .* var.unit
    return Raster(converted, dims(raw2d); crs = crs(raw2d))
end

# ---------------------------------------------------------------------------
# Derivation registry
# ---------------------------------------------------------------------------

# Previously
# const _ENVELOPE_PHYSICS = (
#     Val(:reference_temperature_max),
#     Val(:reference_temperature_min),
#     Val(:mean_temperature),
#     Val(:wind_speed),
#     Val(:actual_vapour_pressure),
#     Val(:vapour_pressure_deficit),
#     Val(:reference_humidity_max),
#     Val(:reference_humidity_min),
#     Val(:reference_wind_max),
#     Val(:reference_wind_min),
#     Val(:cloud_cover),
#     Val(:cloud_min),
#     Val(:cloud_max),
#     Val(:deep_soil_temperature),
# )

const _ENVELOPE_PHYSICS = (
    Temperature(Maximum()),
    Temperature(Minimum()),
    Reference(Temperature(Maximum())),
    Reference(Temperature(Minimum())),
    Temperature(Mean()),
    WindSpeed(),
    ActualVapourPressure(),
    VapourPressureDeficit(),
    Reference(RelativeHumidity(Maximum())),
    Reference(RelativeHumidity(Minimum())),
    Reference(WindSpeed(Maximum())),
    Reference(WindSpeed(Minimum())),
    CloudCover(),
    CloudCover(Minimum()),
    CloudCover(Maximum()),
    SoilTemperature(Mean()),
)

# Previously
# const _NATIVE_PHYSICS = (
#     Val(:solar_geometry),
#     Val(:wind_speed),
#     Val(:actual_vapour_pressure),
# )

const _NATIVE_PHYSICS = (
   WindSpeed(),
   ActualVapourPressure(),
   Pressure(),
)

const _OUTPUT_PHYSICS = (
    Reference(WindSpeed()),
    Reference(RelativeHumidity()),
    SoilTemperature(Mean()),
    Rainfall(DailyTotal()),
)

abstract type ResamplingRule end
struct Interpolate <: ResamplingRule end
struct Disaggregate <: ResamplingRule end
struct Accumulate <: ResamplingRule end

# TODO these dont seem right - how is there enough information just from the kind to know
# the required resampling
@inline resampling_rule(::Sample) = Interpolate()
@inline resampling_rule(::GlobalRadiation) = Disaggregate()
@inline resampling_rule(::Rainfall) = Accumulate()
@inline resampling_rule(r::Reference) = resampling_rule(r.sample)

@inline function _native_at(src, d, i, nin, ndays)
    if i < 0
        d > 1 ? src[(d - 2) * nin + nin] : src[(d - 1) * nin + 1]
    elseif i >= nin
        d < ndays ? src[d * nin + 1] : src[(d - 1) * nin + nin]
    else
        src[(d - 1) * nin + i + 1]
    end
end

function _resample!(out, src, ::Interpolate, nin::Int, nout::Int)
    ndays = length(src) ÷ nin
    @inbounds for d in 1:ndays
        ob = (d - 1) * nout
        for j in 0:(nout - 1)
            x  = ((j + 0.5) / nout) * nin - 0.5   # output midpoint on native index axis
            i0 = floor(Int, x)
            t  = x - i0
            v0 = _native_at(src, d, i0,     nin, ndays)
            v1 = _native_at(src, d, i0 + 1, nin, ndays)
            out[ob + j + 1] = v0 * (1 - t) + v1 * t
        end
    end
    return out
end
function _resample!(out, src, ::Disaggregate, nin::Int, nout::Int, clear_sky)
    ndays = length(src) ÷ nin
    hpb = nout ÷ nin                       # output steps per native block
    @inbounds for d in 1:ndays
        hb = (d - 1) * nout
        for b in 0:(nin - 1)
            cs_sum = zero(eltype(clear_sky))
            for j in 1:hpb
                cs_sum += clear_sky[hb + b * hpb + j]
            end
            cs_mean = cs_sum / hpb
            obs = src[(d - 1) * nin + b + 1]
            opacity = cs_mean <= zero(cs_mean) ? 0.0 :
                clamp(ustrip(u"W/m^2", obs) / ustrip(u"W/m^2", cs_mean), 0.0, 1.0)
            for j in 1:hpb
                out[hb + b * hpb + j] = opacity * clear_sky[hb + b * hpb + j]
            end
        end
    end
    return out
end
function _resample!(daily, src, ::Accumulate, nin::Int)
    ndays = length(daily)
    @inbounds for d in 1:ndays
        s = zero(eltype(src))
        for i in 1:nin
            s += src[(d - 1) * nin + i]
        end
        daily[d] = s
    end
    return daily
end

@inline function _resample_to_output!(buffers, nin::Int, nout::Int, ctx)
    native = buffers.native
    map(_SERIES_OUTPUT) do q
        _resample_quantity!(buffers, native, q, nin, nout, ctx)
    end
    return nothing
end

@inline function _resample_quantity!(buffers, native, q::Sample, nin, nout, ctx)
    name = canonical_name(q)
    hasproperty(native, name) || return nothing
    _resample_by_rule!(getproperty(buffers, name), getproperty(native, name),
                       resampling_rule(q), nin, nout, ctx)
    return nothing
end

@inline function _resample_by_rule!(out, src, rule::Union{Interpolate,Disaggregate},
                                    nin, nout, ctx)
    nin == nout ? copyto!(out, src) : _disaggregate!(out, src, rule, nin, nout, ctx)
    return nothing
end
@inline _disaggregate!(out, src, rule::Interpolate, nin, nout, ctx) =
    _resample!(out, src, rule, nin, nout)
@inline _disaggregate!(out, src, rule::Disaggregate, nin, nout, ctx) =
    _resample!(out, src, rule, nin, nout, ctx.scratch.solar.out.global_horizontal)
@inline function _resample_by_rule!(out, src, ::Accumulate, nin, nout, ctx)
    # Coarsening accumulation (nin != nout) isn't implemented yet -- only the
    # same-resolution passthrough, which is what SubDaily->SubDaily needs.
    nin == nout && copyto!(out, src)
    return nothing
end


@inline function _wind_speed!(out, u, v)
    @. out = sqrt(u^2 + v^2)
    return out
end

function _vapour_pressure_from_specific_humidity!(out, q, p)
    @inbounds for k in eachindex(out)
        out[k] = q[k] * p[k] / (0.622 + 0.378 * q[k])
    end
    return out
end

# Per-year mean of `series`, broadcast across each year's `out` entries.
# Year boundaries come from `days_of_year` (wherever it drops instead of
# increasing), so leap years and monthly-representative-day sequences are
# both handled correctly -- unlike splitting by a fixed `length ÷ nyears`,
# which silently drops the remainder whenever years aren't equal length.
function _deep_soil!(out, series, days_of_year)
    n = length(out)
    d = 1
    @inbounds while d <= n
        d_end = d
        while d_end < n && days_of_year[d_end + 1] > days_of_year[d_end]
            d_end += 1
        end
        s = zero(eltype(series))
        for k in d:d_end
            s += series[k]
        end
        m = s / (d_end - d + 1)
        for k in d:d_end
            out[k] = m
        end
        d = d_end + 1
    end
    return out
end

# ---------------------------------------------------------------------------
# Buffer allocation
# ---------------------------------------------------------------------------

"""
    allocate_weather_buffers(source, target, years; days_of_year = _days_of_year(weather_calendar(source), years))

Allocate per-thread weather scratch buffers sized to `days_of_year`. Pass the
actual requested day-of-year vector for sub-yearly runs; defaults to the
full `years` span otherwise.
"""
function allocate_weather_buffers(source::Type, target::Timestep, years;
                                  days_of_year::AbstractVector{Int} = _days_of_year(weather_calendar(source), years))
    cal = weather_calendar(source)
    native = native_timestep(source)
    _allocate_weather_buffers(cal, native, target, source, days_of_year)
end

# Which env timeseries a Sample's data feeds — `EnvHourly` sends the read
# through `buffers.native.<name>`, later resampled to output and sourced
# by `env_hourly`; `EnvDaily` writes straight into `buffers.<name>`, the
# slot `env_daily` reads. Dispatched on the Sample type so the partition
# below resolves through method dispatch, not a runtime branch.
abstract type EnvSlot end
struct EnvHourly <: EnvSlot end
struct EnvDaily  <: EnvSlot end
struct EnvStatic <: EnvSlot end   # no Ti axis -- read separately via weather_grid_elevation, never per-timestep

@inline env_slot(::Sample)                = EnvHourly()
@inline env_slot(::Rainfall{DailyTotal})  = EnvDaily()
@inline env_slot(::SoilTemperature{Mean}) = EnvDaily()
@inline env_slot(::SoilMoisture)          = EnvDaily()
@inline env_slot(::Elevation)             = EnvStatic()
@inline env_slot(v::Variable)             = env_slot(quantity(v))

# Assembly-time parent buffer for a slot.
@inline slot_target(buffers, ::EnvHourly) = buffers.native
@inline slot_target(buffers, ::EnvDaily)  = buffers

# Compile-time partition — keep vs skip decided by method dispatch on the
# (var_slot, wanted_slot) type pair, no boolean comparison.
@inline vars_for(vars::Tuple, slot::EnvSlot) = _vars_for((), vars, slot)
@inline _vars_for(acc::Tuple, ::Tuple{}, ::EnvSlot) = acc
@inline function _vars_for(acc::Tuple, vars::Tuple, slot::EnvSlot)
    v = first(vars)
    return _vars_for_step(acc, v, env_slot(v), slot, Base.tail(vars))
end
@inline _vars_for_step(acc::Tuple, v, ::S, slot::S, rest::Tuple) where {S<:EnvSlot} =
    _vars_for((acc..., v), rest, slot)
@inline _vars_for_step(acc::Tuple, _v, ::EnvSlot, slot::EnvSlot, rest::Tuple) =
    _vars_for(acc, rest, slot)

# All vars except EnvStatic ones (e.g. Elevation), for the MinMax native path.
@inline non_static_vars(vars::Tuple) = _non_static_vars((), vars)
@inline _non_static_vars(acc::Tuple, ::Tuple{}) = acc
@inline function _non_static_vars(acc::Tuple, vars::Tuple)
    v = first(vars)
    return _non_static_vars_step(acc, v, env_slot(v), Base.tail(vars))
end
@inline _non_static_vars_step(acc::Tuple, _v, ::EnvStatic, rest::Tuple) =
    _non_static_vars(acc, rest)
@inline _non_static_vars_step(acc::Tuple, v, ::EnvSlot, rest::Tuple) =
    _non_static_vars((acc..., v), rest)

# Native quantities only — `EnvDaily` ones live in `buffers.<name>`, not here.
# `_full_variables`, not `variables`: buffers need a slot for fallback fields too.
@inline _native_quantities(source::Type) = _merge_unique(
    map(quantity, vars_for(_full_variables(source), EnvHourly())),
    (WindSpeed(), ActualVapourPressure(), Pressure()),
)

# Append each quantity in `extras` to `base` only if `base` doesn't already
# carry a quantity with the same canonical name. Keeps buffer allocation
# collision-free when a source natively declares one of the always-derivable
# quantities (e.g. BARRA's `WindSpeed`).
@inline _merge_unique(base::Tuple, ::Tuple{}) = base
@inline function _merge_unique(base::Tuple, extras::Tuple)
    first_extra = first(extras)
    name = canonical_name(first_extra)
    already = any(x -> canonical_name(x) === name, base)
    new_base = already ? base : (base..., first_extra)
    return _merge_unique(new_base, Base.tail(extras))
end

@inline _zeros(unit, n) = zeros(typeof(0.0 * unit), n)
@inline _zeros(::Nothing, n) = zeros(Float64, n)

@inline _zero_buffers(quantities::Tuple, n) =
    NamedTuple{map(canonical_name, quantities)}(map(q -> _zeros(working_unit(q), n), quantities))

function _environment_daily(n, rainfall, deep_soil_temperature, shade::Vector{Float64})
    DailyTimeseries(;
        shade,
        soil_wetness = zeros(Float64, n),
        surface_emissivity = fill(0.95, n),
        cloud_emissivity = fill(0.95, n),
        rainfall, deep_soil_temperature,
        leaf_area_index = fill(0.1, n),
    )
end

_environment_hourly(out) = HourlyTimeseries(;
    out.pressure, out.reference_temperature, out.reference_humidity,
    reference_wind_speed = out.wind_speed,
    out.global_radiation, out.longwave_radiation, out.cloud_cover,
    out.rainfall, out.zenith_angle,
)

_environment_hourly_stub(n) = HourlyTimeseries(;
    pressure = Fill(atmospheric_pressure(0.0u"m"), n * 24),
    reference_temperature = nothing, reference_humidity = nothing,
    reference_wind_speed = nothing, global_radiation = nothing,
    longwave_radiation = nothing, cloud_cover = nothing,
    rainfall = nothing, zenith_angle = nothing,
)

const _ENVELOPE_BUFFERS = (
    Temperature(Maximum()), Temperature(Minimum()), Temperature(Mean()),
    Reference(Temperature(Minimum())), Reference(Temperature(Maximum())),
    WindSpeed(), EastwardWindSpeed(), NorthwardWindSpeed(),
    Reference(WindSpeed(Minimum())), Reference(WindSpeed(Maximum())),
    VapourPressureDeficit(),
    ActualVapourPressure(),
    SpecificHumidity(), Pressure(),
    Reference(RelativeHumidity(Minimum())), Reference(RelativeHumidity(Maximum())),
    GlobalRadiation(), CloudCover(), CloudCover(Minimum()), CloudCover(Maximum()),
    Rainfall(), SoilTemperature(Mean()), SoilMoisture(),
)

const _MINMAX_FORCING_INPUTS = (
    :reference_temperature_min, :reference_temperature_max,
    :reference_wind_speed_min, :reference_wind_speed_max,
    :reference_humidity_min, :reference_humidity_max,
    :cloud_cover_min, :cloud_cover_max,
)

function _envelope_forcings(variables::Tuple, b)
    valuesource = name -> getproperty(b, name)
    # minmax_forcings derives each _mean from its _min/_max; the field names
    # it takes match this package's own canonical naming exactly.
    standard = Microclimate.minmax_forcings(; NamedTuple{_MINMAX_FORCING_INPUTS}(b)...)
    timed = _timed_forcings(variables, valuesource)
    isempty(timed) && return standard
    base = (; standard.reference_temperature, standard.reference_wind_speed, standard.cloud_cover)
    humidity = haskey(timed, :actual_vapour_pressure) ?
        (; reference_humidity = Derived(
            RelativeHumidityFromVapourPressureAndTemperature(GoffGratch()),
            (:actual_vapour_pressure, :reference_temperature))) :
        (; standard.reference_humidity)
    return merge(base, timed, humidity)
end

function _timed_forcings(variables::Tuple, valuesource)
    Sample = Tuple{TimeOfDay, Symbol}
    groups = Pair{Symbol, Vector{Sample}}[]
    for v in variables
        q = quantity(v)
        t = qualifier(q)
        t isa TimeOfDay || continue
        key = physical_quantity(q)
        i = findfirst(p -> first(p) === key, groups)
        sample = (t, canonical_name(v))
        i === nothing ? push!(groups, key => Sample[sample]) : push!(last(groups[i]), sample)
    end
    isempty(groups) && return (;)
    names = Tuple(first(g) for g in groups)
    forcings = map(groups) do (_, samples)
        sort!(samples; by = s -> Microclimate.nominal_hour(first(s)))
        times = Tuple(first(s) for s in samples)
        curve = _clock_curve(times)
        Microclimate.DielForcing(curve, Tuple(valuesource(last(s)) for s in samples))
    end
    return NamedTuple{names}(Tuple(forcings))
end

function _clock_curve(times::Tuple)
    n = length(times)
    n >= 2 || error("a timed forcing needs at least two samples, got $n")
    shapes = ntuple(i -> Microclimate.Linear(times[i], times[mod1(i + 1, n)]), n)
    return Microclimate.DielCurve(shapes, times)
end

@inline function _variable_buffers(variables::Tuple, n::Int)
    names = map(canonical_name, variables)
    return NamedTuple{names}(map(v -> _zeros(working_unit(quantity(v)), n), variables))
end

const _SERIES_OUTPUT = (
    Reference(Temperature()), Reference(RelativeHumidity()), WindSpeed(), Pressure(),
    CloudCover(), GlobalRadiation(), LongwaveRadiation(), Rainfall(),
    ZenithAngle(), ActualVapourPressure(),
)
const _SERIES_DAILY = (Rainfall(DailyTotal()), SoilTemperature(Mean()), SoilMoisture())

function _allocate_weather_buffers(calendar::Calendar, ::MinMax, ::Timestep, source::Type,
                                   days_of_year::AbstractVector{Int})
    nsteps = length(days_of_year)
    vars = _full_variables(source)
    b = merge(_zero_buffers(_ENVELOPE_BUFFERS, nsteps), _variable_buffers(vars, nsteps))

    forcings = _envelope_forcings(vars, b)
    environment_minmax = _minmax_env_type(calendar)(; forcings)
    shade = zeros(Float64, nsteps)
    environment_daily = _environment_daily(nsteps, b.rainfall, b.deep_soil_temperature, shade)
    environment_hourly = _environment_hourly_stub(nsteps)

    return (; b..., shade, days_of_year,
        environment_minmax, environment_daily, environment_hourly)
end

function _allocate_weather_buffers(::Calendar, native_step::SubDaily, target_step::SubDaily,
                                   source::Type, days_of_year::AbstractVector{Int})
    native_names = _native_quantities(source)
    ndays = length(days_of_year)                  # 365 per year (no-leap calendar)
    nnat  = ndays * samples_per_day(native_step)  # native steps
    nout  = ndays * samples_per_day(target_step)  # output steps

    output = _zero_buffers(_SERIES_OUTPUT, nout)
    daily  = _zero_buffers(_SERIES_DAILY, ndays)
    native = _zero_buffers(native_names, nnat)

    shade = zeros(Float64, ndays)
    environment_daily  = _environment_daily(ndays, daily.rainfall_daily, daily.deep_soil_temperature, shade)
    environment_hourly = _environment_hourly(output)

    return (; output..., native, daily..., shade, days_of_year,
        environment_minmax = nothing, environment_daily, environment_hourly)
end

function assemble_weather!(
    scratch, weather,
    source::Type,
    site, I::Tuple;
    grid_elevation = site.elevation,
    vapour_pressure_method = GoffGratch(),
    lapse_rate_model::LapseRate = EnvironmentalLapseRate(),
    canonical_overrides::NamedTuple = (;),
    wind_reference_height = 2.0u"m",
    target_timestep::Timestep = Hourly(),
)
    buffers = scratch.weather
    atm_pressure = atmospheric_pressure(site.elevation)
    vars = _full_variables(source)
    cal = weather_calendar(source)
    native = native_timestep(source)

    ctx = (;
        site, grid_elevation, lapse_rate_model, vapour_pressure_method,
        atmospheric_pressure = atm_pressure, scratch,
        days_of_year = buffers.days_of_year,
        # Output day count per year, for the deep-soil annual mean.
        steps_per_year = days_per_year(cal),
        native_timestep = native, target_timestep,
        wind_reference_height,
    )

    _assemble!(native, target_timestep, buffers, weather, source, ctx, vars,
               canonical_overrides, I)
    _apply_canonical_overrides!(buffers, canonical_overrides, I)

    environment_hourly = _finalize_environment_hourly(buffers, vars, atm_pressure)

    return (; buffers.environment_minmax, buffers.environment_daily, environment_hourly)
end

@inline function _assemble!(::MinMax, ::Timestep, buffers, weather, source, ctx,
                            variables, overrides, I)
    _read_native!(buffers, weather, non_static_vars(variables), I)
    _run_physics!(buffers, ctx, _ENVELOPE_PHYSICS, variables, overrides)
    return nothing
end
function _assemble!(native_step::SubDaily, target_step::SubDaily, buffers, weather, source,
                    ctx, variables, overrides, I)
    # Per-slot read: `EnvHourly` vars land in `buffers.native.<name>` (later
    # resampled to output); `EnvDaily` declarations (e.g. BARRA's
    # `Rainfall(DailyTotal())` for `:pr`) go straight into `buffers.<name>`,
    # the slot `env_daily` sources.
    _read_native!(slot_target(buffers, EnvHourly()), weather,
                  vars_for(variables, EnvHourly()), I)
    _read_native!(slot_target(buffers, EnvDaily()),  weather,
                  vars_for(variables, EnvDaily()),  I)
    _run_physics!(buffers.native, ctx, _NATIVE_PHYSICS, variables, overrides)
    # The clear-sky solar field is an input the shortwave disaggregation reads,
    # not a derived quantity — compute it before resampling.
    _compute_solar_field!(ctx)
    _resample_to_output!(buffers, samples_per_day(native_step), samples_per_day(target_step), ctx)
    _run_physics!(buffers, ctx, _OUTPUT_PHYSICS, variables, overrides)
    return nothing
end

# TODO: this is a hack just to fill pressure.
@inline function _finalize_environment_hourly(buffers, variables, atm_pressure)
    buffers.environment_minmax === nothing && return buffers.environment_hourly
    _has_native(Val(:pressure), variables)  && return buffers.environment_hourly
    pressure = Fill(atm_pressure, length(buffers.environment_hourly.pressure))
    return setproperties(buffers.environment_hourly, (; pressure))
end

@inline function _read_native!(buffers, weather, variables::Tuple, I)
    unrolled_map(variables) do var
        _read_one_variable!(buffers, weather, var, I)
        nothing
    end
    return nothing
end

function _read_one_variable!(buffers, weather, var::Variable{Name}, I) where {Name}
    target = getproperty(buffers, Name)
    layer = getproperty(weather, Name)
    transform = var.transform
    unit = var.unit
    # Splat the spatial dim tuple from `DimIndices`; `Ti(k)` selects the
    # k-th timestep regardless of how `layer` stores its dims.
    @inbounds for k in eachindex(target)
        target[k] = transform(layer[I..., Ti(k)]) * unit
    end
    return nothing
end

function _run_physics!(buffers, ctx, chain::Tuple, variables::Tuple, overrides::NamedTuple = (;))
    unrolled_map(chain) do v
        n = canonical_name(v)
        is_available = _has_native(Val(n), variables) || _has_override(Val(n), overrides)
        is_available || derive!(v, buffers, ctx)
        nothing
    end
    return nothing
end

@inline function _apply_canonical_overrides!(buffers, overrides::NamedTuple{K}, I) where K
    unrolled_map(K) do name
        _apply_one_override!(buffers, name, getproperty(overrides, name), I)
        nothing
    end
    return nothing
end
@inline _apply_canonical_overrides!(_, ::NamedTuple{()}, _) = nothing

@inline function _apply_one_override!(buffers, name::Symbol, raster, I)
    target = getproperty(buffers, name)
    _copy_override!(target, raster, I)
    return nothing
end

@inline function _copy_override!(target::AbstractVector, raster, I)
    @inbounds if hasdim(raster, Ti)
        for k in eachindex(target)
            target[k] = raster[I..., Ti(k)]
        end
    else
        value = raster[I...]
        for k in eachindex(target)
            target[k] = value
        end
    end
    return nothing
end

@inline _has_override(::Val{N}, ::NamedTuple{K}) where {N, K} = N in K

@inline _has_native(::Val, ::Tuple{}) = false
@inline _has_native(v::Val{N}, vars::Tuple) where {N} =
    _name_matches(v, first(vars)) || _has_native(v, Base.tail(vars))

@inline _name_matches(::Val{N}, ::Variable{N}) where {N} = true
@inline _name_matches(::Val, ::Variable) = false

# ---------------------------------------------------------------------------
# Derivations
# ---------------------------------------------------------------------------

# function derive!(r::Reference, buffers, ctx)
#     derive!(r, r.sample, buffers, ctx)
# end
function derive!(r::Reference{<:Temperature}, buffers, ctx)
    _lapse_correct!(buffers[canonical_name(r)], buffers[canonical_name(r.sample)], ctx)
end
function derive!(s::Temperature{Maximum}, buffers, ctx)
    buffers[canonical_name(s)] .= buffers.mean_temperature .+ buffers.diurnal_temperature_range ./ 2
    return nothing
end
function derive!(s::Temperature{Minimum}, buffers, ctx)
    buffers[canonical_name(s)] .= buffers.mean_temperature .- buffers.diurnal_temperature_range ./ 2
    return nothing
end
function derive!(s::Temperature{Mean}, buffers, ctx)
    buffers[canonical_name(s)] .= (buffers.reference_temperature_max .+ buffers.reference_temperature_min) ./ 2
    return nothing
end
function derive!(s::ActualVapourPressure, buffers, ctx)
    method = ctx.vapour_pressure_method
    avp = buffers[canonical_name(s)]
    if hasproperty(buffers, :dewpoint_temperature)
        @inbounds for k in eachindex(avp)
            avp[k] = vapour_pressure(method, buffers.dewpoint_temperature[k])
        end
    elseif hasproperty(buffers, :humidity)
        @inbounds for k in eachindex(avp)
            avp[k] = buffers.humidity[k] * vapour_pressure(method, buffers.mean_temperature[k])
        end
    elseif hasproperty(buffers, :reference_humidity)
        # Sub-daily-native sources (e.g. BARRA): hourly RH paired with that hour's reference_temperature.
        @inbounds for k in eachindex(avp)
            avp[k] = buffers.reference_humidity[k] * vapour_pressure(method, buffers.reference_temperature[k])
        end
    else
        _vapour_pressure_from_specific_humidity!(avp,
            buffers.specific_humidity, buffers.pressure)
    end
    return nothing
end
function derive!(::VapourPressureDeficit, buffers, ctx)
    method = ctx.vapour_pressure_method
    @inbounds for k in eachindex(buffers.vapour_pressure_deficit)
        saturation = vapour_pressure(method, buffers.mean_temperature[k])
        buffers.vapour_pressure_deficit[k] = saturation - buffers.actual_vapour_pressure[k]
    end
    return nothing
end
function derive!(::Reference{RelativeHumidity{Maximum}}, buffers, ctx)
    _relative_humidity_from_vpd!(buffers.reference_humidity_max,
        buffers.vapour_pressure_deficit, buffers.mean_temperature,
        buffers.reference_temperature_min, ctx.vapour_pressure_method)
    return nothing
end
function derive!(::Reference{RelativeHumidity{Minimum}}, buffers, ctx)
    _relative_humidity_from_vpd!(buffers.reference_humidity_min,
        buffers.vapour_pressure_deficit, buffers.mean_temperature,
        buffers.reference_temperature_max, ctx.vapour_pressure_method)
    return nothing
end
function derive!(::WindSpeed, buffers, ctx)
    _wind_speed!(buffers.wind_speed, buffers.eastward_wind, buffers.northward_wind)
    return nothing
end
# Pressure at the simulated site. When the source declares mean sea level
# pressure (e.g. BARRA's `:psl`), apply the barometric formula at the DEM
# site elevation using each hour's MSLP as the reference — BARRA's grid
# elevation differs from the DEM-derived site one, so the standard-atmosphere
# reference is not appropriate here. When no MSLP is present, fall back to
# the standard atmosphere at site elevation (matches what
# `_finalize_environment_hourly` did inline for MinMax sources).
function derive!(::Pressure, buffers, ctx)
    elevation = ctx.site.elevation
    if hasproperty(buffers, :sea_level_pressure)
        @inbounds for k in eachindex(buffers.pressure)
            buffers.pressure[k] =
                atmospheric_pressure(elevation; reference_pressure = buffers.sea_level_pressure[k])
        end
    else
        p = atmospheric_pressure(elevation)
        @inbounds for k in eachindex(buffers.pressure)
            buffers.pressure[k] = p
        end
    end
    return nothing
end
function derive!(::Reference{WindSpeed{Maximum}}, buffers, ctx)
    shear = _wind_height_correction(ctx.wind_reference_height)
    @. buffers.reference_wind_speed_max = buffers.wind_speed * shear
    return nothing
end
function derive!(::Reference{WindSpeed{Minimum}}, buffers, ctx)
    @. buffers.reference_wind_speed_min = buffers.reference_wind_speed_max * 0.1
    return nothing
end
function derive!(::CloudCover, buffers, ctx)
    (; scratch, site) = ctx
    # `flat_terrain_template` is constant except for elevation/pressure/lat/lon;
    # swap those four in via `setproperties` (zero-cost).
    flat_terrain = setproperties(scratch.cloud_constants.flat_terrain_template,
        (; site.elevation,
           atmospheric_pressure = ctx.atmospheric_pressure,
           site.latitude, site.longitude),
    )
    cloud_from_solar_radiation!(buffers.cloud_cover, scratch.solar.out,
        buffers.global_radiation, flat_terrain, buffers.days_of_year,
        scratch.solar.buffers;
        solar_model = scratch.cloud_constants.solar_model,
        hours = scratch.cloud_constants.hours,
    )
    return nothing
end
function derive!(::CloudCover{Minimum}, buffers, ctx)
    @. buffers.cloud_cover_min = clamp(buffers.cloud_cover * 0.5, 0.0, 1.0)
    return nothing
end
function derive!(::CloudCover{Maximum}, buffers, ctx)
    @. buffers.cloud_cover_max = clamp(buffers.cloud_cover * 2.0, 0.0, 1.0)
    return nothing
end
function derive!(::SoilTemperature{Mean}, buffers, ctx)
    _deep_soil!(buffers.deep_soil_temperature, _air_temperature(buffers), ctx.days_of_year)
    return nothing
end
function derive!(::Rainfall{DailyTotal}, buffers, ctx)
    _resample!(buffers.rainfall_daily, buffers.native.rainfall,
               Accumulate(), samples_per_day(ctx.native_timestep))
    return nothing
end
function derive!(::Reference{WindSpeed{Nothing}}, buffers, ctx)
    shear = _wind_height_correction(ctx.wind_reference_height)
    @. buffers.wind_speed *= shear
    return nothing
end
function derive!(::Reference{RelativeHumidity{Nothing}}, buffers, ctx)
    method = ctx.vapour_pressure_method
    @inbounds for k in eachindex(buffers.reference_humidity)
        saturation = vapour_pressure(method, buffers.reference_temperature[k])
        if saturation <= zero(saturation)
            buffers.reference_humidity[k] = 1.0
        else
            buffers.reference_humidity[k] =
                clamp(buffers.actual_vapour_pressure[k] / saturation, 0.0, 1.0)
        end
    end
    return nothing
end
# Populate the flat-terrain clear-sky solar field in `scratch.solar.out`. This is
# the baseline the shortwave disaggregation divides observed radiation against; it
# is an input to the chain, produced by the solar model, not a weather quantity.
function _compute_solar_field!(ctx)
    (; scratch, site) = ctx
    flat_terrain = setproperties(scratch.cloud_constants.flat_terrain_template,
        (; site.elevation,
           atmospheric_pressure = ctx.atmospheric_pressure,
           site.latitude, site.longitude))
    solar_radiation!(scratch.solar.out, scratch.solar.buffers,
        scratch.cloud_constants.solar_model;
        solar_terrain = flat_terrain,
        days  = ctx.days_of_year,
        hours = scratch.cloud_constants.hours)
    return nothing
end

# ---------------------------------------------------------------------------
# Helpers reused by derivations
# ---------------------------------------------------------------------------

# Power-law wind height correction
_wind_height_correction(z_ref, z_src = 10.0u"m", α = 0.15) = (z_ref / z_src)^α

@inline _air_temperature(b) =
    hasproperty(b, :reference_temperature) ? b.reference_temperature : b.mean_temperature

@inline function _lapse_correct!(out, src, ctx)
    (; site, grid_elevation, lapse_rate_model) = ctx
    Δz = site.elevation - grid_elevation
    if iszero(Δz)
        @. out = src
    else
        lr = lapse_rate(lapse_rate_model)   # scalar K/m — extract before broadcast
        @. out = src - lr * Δz
    end
    return nothing
end

function _relative_humidity_from_vpd!(
    out::AbstractVector,
    vapour_pressure_deficit::AbstractVector,
    mean_temperature::AbstractVector,
    reference_temperature::AbstractVector,
    method,
)
    @inbounds for k in eachindex(out, vapour_pressure_deficit,
                                 mean_temperature, reference_temperature)
        saturation_at_mean = vapour_pressure(method, mean_temperature[k])
        actual = saturation_at_mean - vapour_pressure_deficit[k]
        if actual <= zero(actual)
            out[k] = 0.0
            continue
        end
        saturation_at_reference = vapour_pressure(method, reference_temperature[k])
        out[k] = saturation_at_reference <= zero(saturation_at_reference) ? 1.0 :
                 clamp(actual / saturation_at_reference, 0.0, 1.0)
    end
    return out
end

# ---------------------------------------------------------------------------
# Date-range helpers
# ---------------------------------------------------------------------------

@inline _minmax_env_type(::Monthly) = MonthlyMinMaxEnvironment
@inline _minmax_env_type(::Daily) = DailyMinMaxEnvironment

@inline _days_of_year(::Monthly, years) = repeat(Microclimate.DEFAULT_DAYS, length(years))
# `init=Int[]` forces a real `vcat` even for single-year `years` — otherwise
# `reduce` short-circuits to a bare unconverted UnitRange, which fails
# SolarRadiation.solar_radiation!'s `Vector{<:Real}` check.
@inline _days_of_year(::Daily, years) =
    reduce(vcat, (1:Dates.daysinyear(y) for y in years); init = Int[])

# Convert an old-style integer years range to a full-calendar Date range.
_years_to_dates(years::AbstractRange{Int}) =
    Date(first(years), 1, 1):Day(1):Date(last(years), 12, 31)

# Derive the integer year span that must be loaded from weather sources to
# cover all dates.
_years_from_dates(d::Date) = year(d):year(d)
_years_from_dates(r::AbstractRange{Date})  = minimum(year, r):maximum(year, r)

# Cumulative day offset at which year `years_v[yi]` begins. Setup-time only,
# so allocation doesn't matter here.
function _cumulative_year_days(years_v::Vector{Int}, yi::Int)
    total = 0
    @inbounds for i in 1:(yi - 1)
        total += Dates.daysinyear(years_v[i])
    end
    return total
end

# Normalise user-supplied dates to a sorted Vector{Date}. Feb 29 is kept —
# leap days are simulated like any other day.
_normalise_dates(d::Date) = [d]
function _normalise_dates(r::AbstractRange{Date})
    v = sort(collect(r))
    isempty(v) && error("No simulation dates supplied.")
    return v
end

function _ti_range_for_dates(::Monthly, ::Timestep, years, dates_vec::Vector{Date})
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = (yi_start - 1) * 12 + month(start_d)
    ti_end   = (yi_end   - 1) * 12 + month(end_d)
    days_doy = Int[]
    d = Date(year(start_d), month(start_d), 1)
    stop = Date(year(end_d), month(end_d), 1)
    while d <= stop
        push!(days_doy, Microclimate.DEFAULT_DAYS[month(d)])
        d += Month(1)
    end
    return ti_start, ti_end, days_doy
end
function _ti_range_for_dates(::Daily, native::Timestep, years, dates_vec::Vector{Date})
    spd = samples_per_day(native)
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = _cumulative_year_days(years_v, yi_start) * spd + (dayofyear(start_d) - 1) * spd + 1
    ti_end   = _cumulative_year_days(years_v, yi_end)   * spd + dayofyear(end_d) * spd
    return ti_start, ti_end, dayofyear.(dates_vec)
end

# One "anchor" Date per solver step — used to carry the year when building
# the output DateTime axis. Monthly: 1st of each selected month. Daily: the
# dates_vec itself (already one per calendar day).
function _step_anchor_dates(::Monthly, dates_vec::Vector{Date})
    seen = Set{Tuple{Int,Int}}()
    result = Date[]
    for d in sort(dates_vec)
        k = (year(d), month(d))
        k in seen && continue
        push!(seen, k)
        push!(result, Date(year(d), month(d), 1))
    end
    return result
end
_step_anchor_dates(::Daily, dates_vec::Vector{Date}) = dates_vec

# Monthly forcing → independent representative days (Fortran "monthly mode");
# Daily forcing → continuous run with state carrying day-to-day.
@inline _time_mode(::Monthly) = Microclimate.NonConsecutiveDayMode()
@inline _time_mode(::Daily) = Microclimate.ConsecutiveDayMode()

