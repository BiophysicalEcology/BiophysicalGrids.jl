"""
Abstract type for lapse rate formulations. Dispatch on concrete subtypes to select
the temperature adjustment used in mesoclimate corrections.
"""
abstract type LapseRate end

"""
    EnvironmentalLapseRate <: LapseRate

Standard environmental (average) lapse rate: 6.5 K per 1000 m.
"""
struct EnvironmentalLapseRate <: LapseRate end

"""
    DryAdiabaticLapseRate <: LapseRate

Dry adiabatic lapse rate: 9.8 K per 1000 m. Used for unsaturated rising air.
"""
struct DryAdiabaticLapseRate <: LapseRate end

"""
    SaturatedAdiabaticLapseRate <: LapseRate

Saturated (moist) adiabatic lapse rate: ~6.0 K per 1000 m.
Approximates average conditions; varies with temperature and pressure.
"""
struct SaturatedAdiabaticLapseRate <: LapseRate end

"""
    CustomLapseRate{R} <: LapseRate

User-supplied lapse rate. Provide `rate` as a Unitful K/m quantity.

# Example
```julia
lr = CustomLapseRate(0.007u"K/m")
```
"""
struct CustomLapseRate{R} <: LapseRate
    rate::R
end

"""
    lapse_rate(::LapseRate) → Quantity{K/m}

Return the lapse rate constant for the given lapse rate type.
"""
lapse_rate(::EnvironmentalLapseRate)      = 0.0065u"K/m"
lapse_rate(::DryAdiabaticLapseRate)       = 0.0098u"K/m"
lapse_rate(::SaturatedAdiabaticLapseRate) = 0.006u"K/m"
lapse_rate(lr::CustomLapseRate)           = lr.rate

"""
    lapse_adjust_temperature(T, Δz, lr::LapseRate)

Adjust temperature `T` for an elevation difference `Δz` using lapse rate `lr`.

`Δz` = target elevation − reference elevation. Positive `Δz` (target is higher)
lowers temperature; negative `Δz` (target is lower) raises it.

Accepts scalar or vector `T`, and any Unitful temperature and length units.

# Examples
```julia
lapse_adjust_temperature(20.0u"°C", 1000.0u"m", EnvironmentalLapseRate())  # ≈ 13.5 °C
lapse_adjust_temperature(20.0u"°C", 1000.0u"m", DryAdiabaticLapseRate())   # ≈ 10.2 °C
```
"""
lapse_adjust_temperature(T, Δz, lr::LapseRate) = T - lapse_rate(lr) * Δz

lapse_adjust_temperature(T::AbstractVector, Δz, lr::LapseRate) =
    lapse_adjust_temperature.(T, Ref(Δz), Ref(lr))

"""
    dewpoint_lapse_adjust(T_dew, Δz)

Adjust dewpoint temperature for an elevation difference `Δz`.
Uses the standard moist adiabatic dewpoint lapse rate of ~0.4 K per 100 m (0.004 K/m).
"""
const DEWPOINT_LAPSE_RATE = 0.004u"K/m"

dewpoint_lapse_adjust(T_dew, Δz) = T_dew - DEWPOINT_LAPSE_RATE * Δz

dewpoint_lapse_adjust(T_dew::AbstractVector, Δz) =
    dewpoint_lapse_adjust.(T_dew, Ref(Δz))
