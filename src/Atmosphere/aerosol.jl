# Aerosol optical depth from the Global Aerosol Data Set (GADS).
#
# Koepke, P., Hess, M., Schult, I., and Shettle, E.P. (1997).
# Global Aerosol Data Set. MPI-Report No. 243, Max-Planck-Institut für
# Meteorologie, Hamburg.
#
# The companion NetCDF file (gads.nc) contains pre-computed spectral aerosol
# optical depths on a 5° global grid for two seasons and eight relative-
# humidity levels.  This module reads that file and interpolates to an
# arbitrary location, humidity, month, and wavelength grid.

"""
    get_aerosol_optical_depth(point, relative_humidity, month;
        gads_file   = joinpath(homedir(), "Spatial_Data", "gads", "gads.nc"),
        interpolate = true,
        wavelengths = SolarRadiation.DEFAULT_WAVELENGTHS,
    ) -> Vector{Float64}

Return the spectral aerosol optical depth vector at the requested location,
relative humidity, and time of year, interpolated to `wavelengths`.

# Arguments
- `point`               : GeoInterface-compatible point geometry (e.g. `Point([lon, lat])`)
- `relative_humidity`   : 0–1 fraction
- `month`               : 1–12; used to blend NH summer (July) and winter
                          (January) climatologies via cosine interpolation

# Keywords
- `gads_file`   : path to `gads.nc`; defaults to
                  `\$RASTERDATASOURCES_PATH/gads/gads.nc` (falls back to
                  `homedir()/gads/gads.nc` if the env var is not set)
- `interpolate` : `true` (default) → bilinear spatial interpolation over the
                  surrounding four 5° grid points; `false` → nearest grid point
- `wavelengths` : target wavelength vector (Unitful, nm); defaults to
                  `SolarRadiation.DEFAULT_WAVELENGTHS` (111 points, 290–4000 nm)

# Returns
`Vector{Float64}` of length `length(wavelengths)`, one optical depth per
wavelength.  Values are ≥ 0 (dimensionless).

# Example
```julia
aod = get_aerosol_optical_depth(Point([6.87, 45.92]), 0.60, 7)   # Chamonix, July, 60 % RH
solar_model = SolarProblem(aerosol_optical_depth = aod)
result = simulate_microclimate(solar_terrain, micro_terrain, soil, weather;
                               solar_model)
```
"""
function get_aerosol_optical_depth(
    point             ,
    relative_humidity :: Real,
    month             :: Int;
    gads_file   :: AbstractString = joinpath(get(ENV, "RASTERDATASOURCES_PATH", homedir()), "gads", "gads.nc"),
    interpolate :: Bool           = true,
    wavelengths                   = SolarRadiation.DEFAULT_WAVELENGTHS,
)
    lon, lat = _lonlat(point)
    isfile(gads_file) || error(
        "GADS file not found: \"$gads_file\"\n" *
        "Set the `gads_file` keyword or place gads.nc at the default path.")

    # ── 1. Load the NetCDF into a DimArray ────────────────────────────────────
    da = NCDatasets.NCDataset(gads_file, "r") do ds
        lons    = Float64.(ds["lon"][:])
        lats    = Float64.(ds["lat"][:])
        rh_raw  = Float64.(ds["relhum"][:])      # stored as 0, 50, 70, … 99
        seasons = Float64.(ds["season"][:])      # 0 = summer, 1 = winter
        wls     = Float64.(ds["wavelength"][:])  # nm
        data    = Float64.(ds["OPTDEPTH"][:, :, :, :, :])  # (lon,lat,rh,season,wl)

        DimArray(data, (
            Dim{:lon}(lons),
            Dim{:lat}(lats),
            Dim{:relative_humidity}(rh_raw ./ 100.0),  # convert to 0–1 fraction
            Dim{:season}(seasons),
            Dim{:wavelength}(wls),
        ))
    end

    # ── 2. Season blend (cosine interpolation, peak summer = July = month 7) ──
    #   season_weight = 1 → pure summer, 0 → pure winter
    season_weight = 0.5 * (1.0 - cos(2π * (month - 1) / 12))
    summer = da[season = At(0.0)]   # (lon, lat, relative_humidity, wavelength)
    winter = da[season = At(1.0)]
    da_sw  = season_weight .* summer .+ (1.0 - season_weight) .* winter

    # ── 3. Relative-humidity interpolation ────────────────────────────────────
    rh_clamped = clamp(Float64(relative_humidity), 0.0, 0.99)
    rh_levels  = dims(da_sw, :relative_humidity).val
    da_rh      = _interp_dim(da_sw, :relative_humidity, rh_clamped, rh_levels)
    # da_rh is now (lon, lat, wavelength)

    # ── 4. Spatial selection / interpolation ──────────────────────────────────
    lon_levels = dims(da_rh, :lon).val
    lat_levels = dims(da_rh, :lat).val

    if interpolate
        da_spatial = _bilinear_latlon(da_rh, Float64(lat), Float64(lon),
                                      lat_levels, lon_levels)
    else
        da_spatial = da_rh[lon = Near(Float64(lon)),
                           lat = Near(Float64(lat))]
    end
    # da_spatial is now a 1-D DimArray indexed by wavelength

    gads_wls  = Float64.(dims(da_spatial, :wavelength).val)   # 25 GADS wavelengths (nm)
    gads_aod  = Float64.(parent(da_spatial))                   # 25 optical depth values

    # ── 5. Wavelength interpolation (GADS 25 pts → target grid) ──────────────
    target_wls_nm = ustrip.(u"nm", wavelengths)
    aod = [max(0.0, _linear_interp1(gads_wls, gads_aod, w)) for w in target_wls_nm]

    return aod
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
Interpolate a DimArray linearly along `dim_name` at value `x`.
`levels` must be sorted ascending.
"""
function _interp_dim(da, dim_name::Symbol, x::Float64, levels::AbstractVector{Float64})
    idx = searchsortedlast(levels, x)
    idx = clamp(idx, 1, length(levels) - 1)
    x0, x1 = levels[idx], levels[idx + 1]
    t = (x - x0) / (x1 - x0)
    s0 = da[Dim{dim_name}(At(x0))]
    s1 = da[Dim{dim_name}(At(x1))]
    return (1.0 - t) .* s0 .+ t .* s1
end

"""
Bilinear interpolation in latitude and longitude.
`da` must have dims `:lon`, `:lat`, and at least one other dim.
"""
function _bilinear_latlon(da, latitude::Float64, longitude::Float64,
                          lat_levels::AbstractVector{Float64},
                          lon_levels::AbstractVector{Float64})
    # Wrap longitude to match GADS grid range (−180 to 175)
    lon = mod(longitude + 180.0, 360.0) - 180.0

    # Find bracketing indices
    i_lat = clamp(searchsortedlast(lat_levels, latitude), 1, length(lat_levels) - 1)
    i_lon = clamp(searchsortedlast(lon_levels, lon),      1, length(lon_levels) - 1)

    lat0, lat1 = lat_levels[i_lat],     lat_levels[i_lat + 1]
    lon0, lon1 = lon_levels[i_lon],     lon_levels[i_lon + 1]

    t_lat = (latitude - lat0) / (lat1 - lat0)
    t_lon = (lon      - lon0) / (lon1 - lon0)

    v00 = da[lon = At(lon0), lat = At(lat0)]
    v10 = da[lon = At(lon1), lat = At(lat0)]
    v01 = da[lon = At(lon0), lat = At(lat1)]
    v11 = da[lon = At(lon1), lat = At(lat1)]

    row0 = (1.0 - t_lon) .* v00 .+ t_lon .* v10   # lat0 row
    row1 = (1.0 - t_lon) .* v01 .+ t_lon .* v11   # lat1 row
    return (1.0 - t_lat) .* row0 .+ t_lat .* row1
end

"""
Piecewise-linear 1-D interpolation of `ys` at `x` given sorted `xs`.
Clamps to the endpoint values outside the range.
"""
function _linear_interp1(xs::AbstractVector{Float64}, ys::AbstractVector{Float64},
                         x::Float64)
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return (1.0 - t) * ys[i] + t * ys[i + 1]
end
