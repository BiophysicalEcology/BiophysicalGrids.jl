# Placeholder: point microclimate simulation with ERA5 hourly forcing
#
# ERA5 provides hourly reanalysis data (~31 km) and will populate
# HourlyTimeseries directly rather than MonthlyMinMaxEnvironment.
#
# Planned API (same two-step pattern):
#
#   weather = get_weather(ERA5, lon, lat; tstart = DateTime(2000,1,1), tend = DateTime(2000,12,31))
#   result  = simulate_microclimate(solar_terrain, micro_terrain, soil_model, weather)
