# Terrain utility functions for biophysical grid analyses.
#
# These functions bridge gridded DEM rasters and the SolarRadiation.jl /
# Geomorphometry.jl APIs.  Horizon angles are computed via
# Geomorphometry.horizon_angle (KernelAbstractions-based sweep algorithm).


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
    load_utm_dem(center_lon, center_lat, extent_lon, extent_lat) -> NamedTuple

Download an SRTM DEM tile, crop to the requested bounding box, and reproject
to the local UTM zone.

# Returns
NamedTuple with fields:
- `utm_dem`       : reprojected `Raster`
- `x_coords_utm`  : easting vector (m)
- `y_coords_utm`  : northing vector (m)
- `nx_utm`        : number of columns
- `ny_utm`        : number of rows
- `cs`            : `(dx, dy)` cell size tuple (m)
"""
function load_utm_dem(center_lon, center_lat, extent_lon, extent_lat)
    lon_min = center_lon - extent_lon / 2
    lon_max = center_lon + extent_lon / 2
    lat_min = center_lat - extent_lat / 2
    lat_max = center_lat + extent_lat / 2

    tile_paths  = getraster(SRTM; bounds = (lon_min, lat_min, lon_max, lat_max))
    valid_paths = filter(!ismissing, vec(tile_paths))
    isempty(valid_paths) && error("No SRTM tile found for bounds " *
        "($(lon_min), $(lat_min), $(lon_max), $(lat_max)).")
    dem_full  = Raster(only(valid_paths))
    dem_wgs84 = dem_full[X(Between(lon_min, lon_max)), Y(Between(lat_min, lat_max))]

    utm_crs      = get_utm_crs(dem_wgs84)
    utm_dem      = Rasters.resample(dem_wgs84; crs = utm_crs, method = :bilinear)
    x_coords_utm = collect(lookup(utm_dem, X))
    y_coords_utm = collect(lookup(utm_dem, Y))
    nx_utm       = length(x_coords_utm)
    ny_utm       = length(y_coords_utm)
    cs           = (abs(x_coords_utm[2] - x_coords_utm[1]),
                    abs(y_coords_utm[2] - y_coords_utm[1]))

    return (; utm_dem, x_coords_utm, y_coords_utm, nx_utm, ny_utm, cs)
end

"""
    compute_terrain_grids(utm_dem, x_coords_utm, y_coords_utm; n_horizon_angles=32)
        -> NamedTuple

Compute the full set of per-pixel terrain grids needed for solar radiation and
microclimate simulations from a UTM-projected DEM raster.

Handles the (X,Y) vs (Y,X) internal Rasters dimension ordering automatically.

# Returns
NamedTuple with fields:
- `dem_data`      : `(ny, nx)` Float64 elevation matrix (NaN for no-data)
- `data_is_xy`    : `true` if the raster stores data as `(nx, ny)`
- `y_descending`  : `true` if northing decreases along the first axis
- `elevation_m`   : elevation raster with `u"m"` units
- `slope_deg`     : slope raster with `u"°"` units
- `aspect_deg`    : aspect raster with `u"°"` units
- `latitude_deg`  : per-pixel latitude raster with `u"°"` units
- `longitude_deg` : per-pixel longitude raster with `u"°"` units
- `pressure_r`    : atmospheric pressure raster (`Pa`)
- `horizons_u`    : `(ny, nx, n_dirs)` horizon-angle array with `u"°"` units
"""
function compute_terrain_grids(utm_dem, x_coords_utm, y_coords_utm;
                               n_horizon_angles = 32)
    nx_utm = length(x_coords_utm)
    ny_utm = length(y_coords_utm)
    cs     = (abs(x_coords_utm[2] - x_coords_utm[1]),
              abs(y_coords_utm[2] - y_coords_utm[1]))

    # ---- Internal dimension order ----------------------------------------
    dem_data_raw = map(x -> ismissing(x) ? NaN : Float64(x), utm_dem.data)
    data_is_xy   = size(dem_data_raw) == (nx_utm, ny_utm) && ny_utm != nx_utm
    dem_data     = data_is_xy ? permutedims(dem_data_raw) : dem_data_raw  # (ny, nx)

    # ---- Slope & aspect (Geomorphometry needs (nx,ny), Y ascending) -------
    dem_for_geomorph = data_is_xy ? dem_data_raw : permutedims(dem_data)
    y_descending     = y_coords_utm[1] > y_coords_utm[end]
    if y_descending
        dem_for_geomorph = reverse(dem_for_geomorph, dims = 2)
    end

    slope_g  = Geomorphometry.slope( dem_for_geomorph; method = Horn(), cellsize = cs)
    aspect_g = Geomorphometry.aspect(dem_for_geomorph; method = Horn(), cellsize = cs)

    slope_out  = y_descending ? reverse(slope_g,  dims = 2) : slope_g
    aspect_out = y_descending ? reverse(aspect_g, dims = 2) : aspect_g

    if data_is_xy
        utm_slope  = Raster(slope_out;  dims = dims(utm_dem))
        utm_aspect = Raster(aspect_out; dims = dims(utm_dem))
    else
        utm_slope  = Raster(permutedims(slope_out);  dims = dims(utm_dem))
        utm_aspect = Raster(permutedims(aspect_out); dims = dims(utm_dem))
    end

    # ---- Per-pixel lat/lon -----------------------------------------------
    latlon_dem            = Rasters.resample(utm_dem; crs = GeoFormatTypes.EPSG(4326))
    lat_south, lat_north  = extrema(lookup(latlon_dem, Y))
    lon_west,  lon_east   = extrema(lookup(latlon_dem, X))

    lats = y_coords_utm[1] < y_coords_utm[end] ?
        range(lat_south, lat_north; length = ny_utm) :
        range(lat_north, lat_south; length = ny_utm)
    lons = x_coords_utm[1] < x_coords_utm[end] ?
        range(lon_west, lon_east; length = nx_utm) :
        range(lon_east, lon_west; length = nx_utm)

    lat_matrix = repeat(collect(lats), 1, nx_utm)
    lon_matrix = repeat(collect(lons)', ny_utm, 1)

    if data_is_xy
        lat_raster = Raster(permutedims(lat_matrix), dims(utm_dem))
        lon_raster = Raster(permutedims(lon_matrix), dims(utm_dem))
    else
        lat_raster = Raster(lat_matrix, dims(utm_dem))
        lon_raster = Raster(lon_matrix, dims(utm_dem))
    end

    # ---- Horizon angles (Geomorphometry.jl) -----------------------------
    # dem_for_geomorph is (nx, ny) with Y ascending — same layout required by
    # slope/aspect above.  horizon_angle returns (nx, ny, n_dirs); undo the Y
    # flip if applied, then permute to (ny, nx, n_dirs) to match the pipeline.
    horizons_g = Geomorphometry.horizon_angle(dem_for_geomorph;
                     directions = n_horizon_angles, cellsize = cs)
    if y_descending
        horizons_g = reverse(horizons_g, dims = 2)
    end
    horizons_u = permutedims(horizons_g, (2, 1, 3)) .* 1.0u"°"

    # ---- Unit-tagged rasters --------------------------------------------
    tag(r, u) = map(x -> (ismissing(x) || isnan(x)) ? missing : x * u, r)

    elevation_m   = tag(utm_dem,    1.0u"m")
    slope_deg     = tag(utm_slope,  1.0u"°")
    aspect_deg    = tag(utm_aspect, 1.0u"°")
    latitude_deg  = tag(lat_raster, 1.0u"°")
    longitude_deg = tag(lon_raster, 1.0u"°")
    pressure_r    = map(e -> ismissing(e) ? missing : atmospheric_pressure(e), elevation_m)

    return (; dem_data, data_is_xy, y_descending,
              elevation_m, slope_deg, aspect_deg,
              latitude_deg, longitude_deg, pressure_r,
              horizons_u)
end

"""
    ascending_y(y_coords, data) -> (y_asc, data_asc)

Return `y_coords` and the `(ny, nx)` matrix `data` both flipped so that y
increases from first to last row, as required by Plots.jl `heatmap`.
"""
function ascending_y(y_coords, data)
    if y_coords[1] > y_coords[end]
        return reverse(y_coords), data[end:-1:1, :]
    else
        return y_coords, data
    end
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

# Returns
`Array{Float64,3}` of size `(ny, nx, n_dirs)` with horizon elevation angles
in degrees.  `NaN` where `dem` is `NaN`.
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
