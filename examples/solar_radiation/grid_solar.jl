# Placeholder: solar radiation simulation across a spatial grid
#
# Will use SolarRadiation.jl to compute clear-sky radiation for every
# grid cell in a raster, accounting for slope, aspect, and horizon angles
# derived from a DEM via Geomorphometry.jl.
#
# Planned API:
#
#   dem     = getraster(SRTM30, ...; extent)
#   results = solar_radiation_grid(dem; year = 2000, hours = 0:23)
