"""
    wind_profile_adjust(ws, z_target, z_ref, roughness_height)

Adjust wind speed `ws` measured at height `z_ref` to height `z_target` using the
logarithmic wind profile (neutral atmospheric stability).

All height arguments must have compatible Unitful units (e.g. metres).

# Example
```julia
# Convert 2 m wind speed to 10 m reference height
ws_10m = wind_profile_adjust(3.0u"m/s", 10.0u"m", 2.0u"m", 0.01u"m")
```
"""
function wind_profile_adjust(ws, z_target, z_ref, roughness_height)
    return ws * log(z_target / roughness_height) / log(z_ref / roughness_height)
end
