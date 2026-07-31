# SLGA -- same 6 depth bins as SoilGrids, already plain physical units.
# No point-query API in PointDataSources.jl for this source; point queries
# below reuse the extent-based path with a small buffer, like SoilGrids.

texture_variables(::Type{SLGA}) = (
    TextureVariable(:bulk_density, :bdod, u"Mg/m^3"),  # already g/cm^3 == Mg/m^3
    TextureVariable(:clay, :clay),
    TextureVariable(:silt, :silt),
    TextureVariable(:sand, :sand),
)

function _load_soil_texture_native(::Type{SLGA}, area::Extent; component = "EV")
    vars = texture_variables(SLGA)
    depth_bins = collect(depths(SLGA))  # getraster's depth::AbstractArray dispatch needs a Vector, not a Tuple
    values = map(vars) do var
        paths = getraster(SLGA, native_field(var); depth = depth_bins, component)
        _texture_values_from_paths(paths, area, var)
    end
    return NamedTuple{map(canonical_name, vars)}(values)
end

function _load_soil_texture_native(::Type{SLGA}, lon::Real, lat::Real; component = "EV")
    area = Extent(
        X = (lon - _SOIL_POINT_BUFFER_DEG, lon + _SOIL_POINT_BUFFER_DEG),
        Y = (lat - _SOIL_POINT_BUFFER_DEG, lat + _SOIL_POINT_BUFFER_DEG),
    )
    return _load_soil_texture_native(SLGA, area; component)
end
