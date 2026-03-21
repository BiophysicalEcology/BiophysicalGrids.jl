# Terrain utility functions for biophysical grid analyses.
#
# These functions bridge gridded DEM rasters and the SolarRadiation.jl /
# Geomorphometry.jl APIs.  The horizon-angle algorithm here is a temporary
# placeholder: once Geomorphometry.jl adds its own horizon routine it will
# be preferred and this function removed.

"""
    get_utm_crs(raster) -> GeoFormatTypes.EPSG

Return the UTM EPSG code whose zone contains the centre of `raster`.
Northern-hemisphere zones are EPSG 326xx; southern-hemisphere zones 327xx.
"""
function get_utm_crs(raster)
    center_lon = lookup(raster, X)[size(raster, X) ÷ 2]
    center_lat = lookup(raster, Y)[size(raster, Y) ÷ 2]
    utm_zone   = floor(Int, (center_lon + 180) / 6) + 1
    epsg_code  = center_lat ≥ 0 ? "326$(lpad(utm_zone, 2, '0'))" :
                                   "327$(lpad(utm_zone, 2, '0'))"
    return GeoFormatTypes.EPSG(parse(Int, epsg_code))
end

"""
    compute_horizon_angles(dem, x_coords, y_coords, n_dirs; verbose=true)
        -> Array{Float64,3}

Compute the terrain horizon elevation angle (degrees above horizontal) for
`n_dirs` evenly-spaced azimuth directions.

# Arguments
- `dem`       : `(ny, nx)` matrix of elevations in consistent linear units
                (usually metres after UTM reprojection).  `NaN` for no-data.
- `x_coords`  : x-axis coordinate vector of length `nx` (easting in metres).
- `y_coords`  : y-axis coordinate vector of length `ny` (northing in metres).
- `n_dirs`    : number of azimuth directions (e.g. 24 for 15° steps).
- `verbose`   : print progress to stdout (default `true`).

# Returns
Array of shape `(ny, nx, n_dirs)` with horizon angles in **degrees**.
Direction `d` corresponds to azimuth `360° × (d-1) / n_dirs` measured
clockwise from north.

# Note
This is a pure-Julia ray-marching implementation intended as a placeholder
until a native Geomorphometry.jl horizon routine becomes available.
"""
function compute_horizon_angles(dem, x_coords, y_coords, n_dirs; verbose = true)
    ny, nx   = size(dem)
    horizons = zeros(Float64, ny, nx, n_dirs)
    dx       = Float64(x_coords[2] - x_coords[1])
    dy       = Float64(y_coords[2] - y_coords[1])

    for d in 1:n_dirs
        az      = 2π * (d - 1) / n_dirs
        raw_dj  = sin(az) / abs(dx)
        raw_di  = cos(az) / dy
        mc      = max(abs(raw_di), abs(raw_dj))
        mc ≈ 0  && continue
        step_di = raw_di / mc
        step_dj = raw_dj / mc

        verbose && print("  horizon direction $d/$n_dirs  \r")
        for j in 1:nx, i in 1:ny
            if isnan(dem[i, j])
                horizons[i, j, d] = NaN
                continue
            end
            elev0   = dem[i, j]
            max_tan = 0.0
            k       = 1
            while true
                ri = round(Int, i + k * step_di)
                rj = round(Int, j + k * step_dj)
                (ri < 1 || ri > ny || rj < 1 || rj > nx) && break
                if !isnan(dem[ri, rj])
                    dist = sqrt((x_coords[rj] - x_coords[j])^2 +
                                (y_coords[ri] - y_coords[i])^2)
                    if dist > 0
                        t = (dem[ri, rj] - elev0) / dist
                        t > max_tan && (max_tan = t)
                    end
                end
                k += 1
            end
            horizons[i, j, d] = atand(max_tan)
        end
    end
    verbose && println()
    return horizons
end
