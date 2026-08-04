"""
    PointQuery

Opt-in `Loader` that fetches each layer via `PointDataSources.getpoint`
(OPeNDAP/REST point queries) instead of downloading/cropping a raster grid.
Not any source's default -- select it per run with e.g.
`MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()`.
Area-mode (`_load_layers`) is single-point only (`area` must be point-like);
points-mode (`_load_layers_at_points`) queries every point in `points_dim`.
"""
struct PointQuery <: Loader end

function _load_layers(::PointQuery, source, fields::Tuple, area::Extent, years)
    width_x, width_y = area.X[2] - area.X[1], area.Y[2] - area.Y[1]
    max(width_x, width_y) < 1e-6 || throw(ArgumentError(
        "PointQuery requires a point-like area (got X width $width_x, Y width $width_y) -- " *
        "it only ever loads one location and would silently mis-load a genuine multi-point run"
    ))
    lon, lat = (area.X[1] + area.X[2]) / 2, (area.Y[1] + area.Y[2]) / 2
    date_range = (Date(first(years), 1, 1), Date(last(years), 12, 31))
    extras = _extra_getpoint_kwargs(source)
    layers = map(fields) do name
        @info "  querying $source $name at ($lon, $lat)..."
        nt = PointDataSources.getpoint(source, name; lon, lat, date=date_range, extras...)
        Raster(reshape(nt.values, 1, 1, length(nt.values)), (X([lon]), Y([lat]), Ti(nt.times)))
    end
    return NamedTuple{fields}(layers)
end

function _load_layers_at_points(::PointQuery, source, fields::Tuple, points_dim, years)
    coords = lookup(points_dim)
    date_range = (Date(first(years), 1, 1), Date(last(years), 12, 31))
    extras = _extra_getpoint_kwargs(source)
    layers = map(fields) do name
        @info "  querying $source $name ($(length(coords)) point(s))..."
        series = [PointDataSources.getpoint(source, name; lon, lat, date=date_range, extras...)
                  for (lon, lat) in coords]
        times = series[1].times
        data = permutedims(reduce(hcat, [s.values for s in series]))
        Raster(data, (points_dim, Ti(times)))
    end
    return NamedTuple{fields}(layers)
end
