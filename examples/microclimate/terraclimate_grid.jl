# Placeholder: grid simulation with TerraClimate
#
# Future implementation will iterate simulate_microclimate over a raster extent
# or a set of grid-cell coordinates, returning a stack of MicroResult objects
# (or a summary raster).
#
# Planned API:
#
#   extent = Extents.Extent(X = (-90.0, -88.0), Y = (42.0, 44.0))
#   results = simulate_microclimate(TerraClimate, extent; ystart = 2000, ...)
#
# Grid runs use Threads.@threads over pixels for parallelism.
