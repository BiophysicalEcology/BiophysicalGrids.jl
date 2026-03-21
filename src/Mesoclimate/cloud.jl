"""
    cloud_from_srad(srad, solar_terrain, doy; solar_model, hours)

Estimate fractional cloud cover (0–1) from observed mean daily solar radiation `srad`
(W/m²) by comparing with clear-sky daily mean radiation from `SolarRadiation.solar_radiation`.

`srad` is the observed mean shortwave radiation over 24 hours (including night = 0),
as provided by gridded datasets such as TerraClimate. The clear-sky baseline is computed
using `solar_model` and `solar_terrain` at day-of-year `doy`.

Returns a value clamped to [0, 1]. Returns 1.0 (full cloud cover) if clear-sky radiation
is zero (polar night or numerical issue).

# Example
```julia
cloud = cloud_from_srad(150.0u"W/m^2", solar_terrain, 105)
```
"""
function cloud_from_srad(
    srad,
    solar_terrain::SolarTerrain,
    doy::Int;
    solar_model::SolarProblem = SolarProblem(),
    hours = collect(0.0:1.0:23.0),
)
    sol = solar_radiation(solar_model; solar_terrain, days = [doy], hours)
    clearsky_mean = sum(sol.global_horizontal) / length(hours)
    clearsky_mean <= zero(clearsky_mean) && return 1.0
    cloud = 1.0 - ustrip(u"W/m^2", srad) / ustrip(u"W/m^2", clearsky_mean)
    return clamp(cloud, 0.0, 1.0)
end

"""
    cloud_from_srad(srad_vec, solar_terrain, doys; kwargs...)

Vector version: estimate cloud cover for multiple day-of-year values.
`srad_vec` and `doys` must be the same length.
"""
function cloud_from_srad(
    srad_vec::AbstractVector,
    solar_terrain::SolarTerrain,
    doys::AbstractVector{Int};
    kwargs...,
)
    return [cloud_from_srad(s, solar_terrain, d; kwargs...) for (s, d) in zip(srad_vec, doys)]
end
