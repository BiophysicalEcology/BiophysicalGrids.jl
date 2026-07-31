# Midpoints of the shared 6 depth bins (0-5,5-15,15-30,30-60,60-100,100-200cm).
const _TEXTURE_DEPTH_MIDPOINTS_CM = [2.5, 10.0, 22.5, 45.0, 80.0, 150.0]

function _interpolate_onto_depths(known_values::AbstractVector{<:Real}, known_depths_cm::AbstractVector{<:Real},
        query_depths_cm::AbstractVector{<:Real})
    spline = CubicSpline(known_values, known_depths_cm; extrapolation = ExtrapolationType.Extension)
    return spline.(query_depths_cm)
end

function _warn_if_texture_incomplete(clay, silt, sand)
    for i in eachindex(clay, silt, sand)
        total = clay[i] + silt[i] + sand[i]
        isapprox(total, 100; atol = 5) ||
            @warn "clay+silt+sand = $total at texture depth $i (expected ~100) -- check units/source."
    end
    return nothing
end

# Matches `Microclimate.example_soil_profile`'s scalar-or-vector idiom.
_broadcast_over_depths(value::AbstractVector, _depths) = value
_broadcast_over_depths(value, depths) = fill(value, length(depths))

function _resolve_root_density(::Nothing, depths)
    default = Microclimate.example_campbell_hydraulic_profile().root_density
    depths == Microclimate.DEFAULT_DEPTHS && return default
    default_depths_cm = ustrip.(u"cm", Microclimate.DEFAULT_DEPTHS)
    query_depths_cm = ustrip.(u"cm", depths)
    return _interpolate_onto_depths(ustrip.(u"m/m^3", default), default_depths_cm, query_depths_cm) .* u"m/m^3"
end
_resolve_root_density(root_density, depths) = _broadcast_over_depths(root_density, depths)

function _assemble_soil_profile(native::NamedTuple;
        depths = Microclimate.DEFAULT_DEPTHS,
        pedotransfer_model::PedotransferModel = CosbyMultivariate(),
        mineral_density = 2.560u"Mg/m^3",
        mineral_conductivity = 1.25u"W/m/K",
        mineral_heat_capacity = 870.0u"J/kg/K",
        root_density = nothing)
    _warn_if_texture_incomplete(native.clay, native.silt, native.sand)

    query_depths_cm = ustrip.(u"cm", depths)
    clay = _interpolate_onto_depths(native.clay, _TEXTURE_DEPTH_MIDPOINTS_CM, query_depths_cm)
    silt = _interpolate_onto_depths(native.silt, _TEXTURE_DEPTH_MIDPOINTS_CM, query_depths_cm)
    sand = _interpolate_onto_depths(native.sand, _TEXTURE_DEPTH_MIDPOINTS_CM, query_depths_cm)
    bulk_density_cm = ustrip.(u"Mg/m^3", native.bulk_density)
    bulk_density = _interpolate_onto_depths(bulk_density_cm, _TEXTURE_DEPTH_MIDPOINTS_CM, query_depths_cm) .* u"Mg/m^3"

    per_depth = map(clay, silt, sand, bulk_density) do clay_i, silt_i, sand_i, bulk_density_i
        (; pedotransfer(pedotransfer_model, clay_i, silt_i, sand_i, bulk_density_i)...,
           field_capacity_and_wilting_point(clay_i, silt_i)...)
    end
    campbell_b = map(p -> p.campbell_b, per_depth)
    air_entry_potential = map(p -> p.air_entry_potential, per_depth)
    saturated_conductivity = map(p -> p.saturated_conductivity, per_depth)
    field_capacity = map(p -> p.field_capacity, per_depth)
    wilting_point = map(p -> p.wilting_point, per_depth)

    soil_profile = SoilProfile(;
        bulk_density,
        mineral_density = _broadcast_over_depths(mineral_density, depths),
        mineral_conductivity = _broadcast_over_depths(mineral_conductivity, depths),
        mineral_heat_capacity = _broadcast_over_depths(mineral_heat_capacity, depths),
        hydraulics = CampbellHydraulicProfile(;
            air_entry_water_potential = air_entry_potential,
            saturated_hydraulic_conductivity = saturated_conductivity,
            campbell_b_parameter = campbell_b,
            root_density = _resolve_root_density(root_density, depths),
        ),
    )
    return (; soil_profile, campbell_b, air_entry_potential, saturated_conductivity, field_capacity, wilting_point)
end

"""
    build_soil_profile(source, area_or_point; depths, pedotransfer_model,
                       mineral_density, mineral_conductivity, mineral_heat_capacity,
                       root_density, quantile, component)
        -> (; soil_profile, campbell_b, air_entry_potential, saturated_conductivity,
              field_capacity, wilting_point)

Fetch soil texture for `source` (`SoilGrids` or `SLGA`) over `area_or_point`
(an `Extent`, spatially averaged into one uniform profile; or `(lon, lat)`
for a point query, buffered into a small `Extent`), interpolate onto `depths`, run
`pedotransfer_model`, and assemble a `SoilProfile`.

`mineral_density`/`mineral_conductivity`/`mineral_heat_capacity`/`root_density`
aren't produced by pedotransfer; `root_density` defaults to
`Microclimate.example_campbell_hydraulic_profile`'s table, re-interpolated
onto `depths` if needed. `campbell_b`/`air_entry_potential`/
`saturated_conductivity` (the per-depth pedotransfer output, already folded
into `soil_profile.hydraulics`) and `field_capacity`/
`wilting_point` are also returned directly, for reference.
"""
function build_soil_profile(::Type{SoilGrids}, area::Extent; quantile = "mean", kw...)
    native = _load_soil_texture_native(SoilGrids, area; quantile)
    return _assemble_soil_profile(native; kw...)
end

function build_soil_profile(::Type{SoilGrids}, point::Tuple{<:Real, <:Real}; quantile = "mean", kw...)
    lon, lat = point
    native = _load_soil_texture_native(SoilGrids, lon, lat; quantile)
    return _assemble_soil_profile(native; kw...)
end

function build_soil_profile(::Type{SLGA}, area::Extent; component = "EV", kw...)
    native = _load_soil_texture_native(SLGA, area; component)
    return _assemble_soil_profile(native; kw...)
end

function build_soil_profile(::Type{SLGA}, point::Tuple{<:Real, <:Real}; component = "EV", kw...)
    lon, lat = point
    native = _load_soil_texture_native(SLGA, lon, lat; component)
    return _assemble_soil_profile(native; kw...)
end
