# Ectotherm thermoregulation at Madison, WI using TerraClimate year 2000.
#
# Demonstrates the full pipeline:
#   1. get_weather / apply_climate_scenario  — TerraClimate forcing (BiophysicalGrids)
#   2. simulate_microclimate × 2             — 0% and 90% shade (BiophysicalGrids)
#   3. thermoregulate                        — ectotherm behaviour (BiophysicalBehaviour)
#
# The organism is a generic 20 g lizard with temperate-zone temperature tolerances.
# Body allometry uses the DesertIguana shape (a cylindrical lizard form).
#
# Requires BiophysicalBehaviour.jl (not a core BiophysicalGrids dependency):
#   using Pkg; Pkg.add(url="https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl")

using BiophysicalGrids
using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Microclimate          # for DailyTimeseries
using SolarRadiation
using FluidProperties
using Unitful
using Statistics
using Plots
using StatsPlots

# ── Location ──────────────────────────────────────────────────────────────
lon, lat  = -89.4557, 43.1379
elevation = 270.0u"m"

months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
hours  = collect(0.0:1.0:23.0)
ndays  = 12
nsteps = ndays * 24

minimum_shade = 0.0
maximum_shade = 0.9

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [1.0, 50.0, 100.0]u"cm"    # ground node + two climbing heights

# ── Step 1: TerraClimate weather ──────────────────────────────────────────
println("Obtaining TerraClimate weather (year 2000)...")
year = 2000
weather          = get_weather(TerraClimate, lon, lat; ystart = year, elevation)
weather_scenario = apply_climate_scenario(TerraClimate{Historical}, ystart = year, weather, lon, lat)

# ── Step 2: Terrain and soil ──────────────────────────────────────────────
solar_terrain = SolarTerrain(;
    elevation,
    slope                = 0.0u"°",
    aspect               = 0.0u"°",
    horizon_angles       = fill(0.0u"°", 24),
    albedo               = 0.15,
    atmospheric_pressure = atmospheric_pressure(elevation),
    latitude             = lat * u"°",
    longitude            = lon * u"°",
)

micro_terrain = MicroTerrain(;
    elevation,
    roughness_height = 0.004u"m",
    karman_constant  = 0.4,
    dyer_constant    = 16.0,
    viewfactor       = 1.0,
)

soil_thermal_model = CampbelldeVriesSoilThermal(;
    bulk_density          = 1.3u"Mg/m^3",
    mineral_density       = 2.56u"Mg/m^3",
    de_vries_shape_factor = 0.1,
    mineral_conductivity  = 2.5u"W/m/K",
    mineral_heat_capacity = 870.0u"J/kg/K",
    saturation_moisture   = 0.42u"m^3/m^3",
    recirculation_power   = 4.0,
    return_flow_threshold = 0.162,
)

# ── Step 3: Microclimate at 0% and 90% shade ──────────────────────────────
# Replace the shade fraction in weather.environment_daily, keeping all other fields.
function with_shade(w, shade_frac)
    ed = w.environment_daily
    n  = length(ed.shade)
    new_ed = DailyTimeseries(;
        shade                 = fill(shade_frac, n),
        soil_wetness          = ed.soil_wetness,
        surface_emissivity    = ed.surface_emissivity,
        cloud_emissivity      = ed.cloud_emissivity,
        rainfall              = ed.rainfall,
        deep_soil_temperature = ed.deep_soil_temperature,
        leaf_area_index       = ed.leaf_area_index,
    )
    merge(w, (; environment_daily = new_ed))
end

common_micro_kwargs = (;
    depths,
    heights,
    runmoist         = false,
    organic_soil_cap = true,
)

println("Solving microclimate (0% shade)...")
low_shade_result = simulate_microclimate(
    solar_terrain, micro_terrain, soil_thermal_model,
    with_shade(weather_scenario, minimum_shade);
    common_micro_kwargs...,
)

println("Solving microclimate (90% shade)...")
high_shade_result = simulate_microclimate(
    solar_terrain, micro_terrain, soil_thermal_model,
    with_shade(weather_scenario, maximum_shade);
    common_micro_kwargs...,
)

available_environments = AvailableEnvironments(
    low_shade_result, high_shade_result, minimum_shade, maximum_shade, depths, heights
)

# ── Step 4: Organism ──────────────────────────────────────────────────────
# Generic 20 g lizard. Temperature tolerances approximate a temperate-zone species.
body = Body(DesertIguana(20.0u"g", 1000.0u"kg/m^3"), Naked())

organism_traits = example_ectotherm_organism_traits(
    activity_period         = Diurnal(),
    T_target                = u"K"(32.0u"°C"),
    T_active_min            = u"K"(28.0u"°C"),
    T_active_max            = u"K"(37.0u"°C"),
    T_bask                  = u"K"(20.0u"°C"),
    T_emerge                = u"K"(15.0u"°C"),
    T_critical_min          = u"K"(5.0u"°C"),
    T_critical_max          = u"K"(42.0u"°C"),
    can_climb               = true,
    can_retreat_underground = true,
    depth_min_underground   = 2,
    burrow_shade_mode       = MaxShadeOnly(),
    warm_signal             = 0.0u"K/hr",
    can_seek_shade          = true,
    shade_min               = minimum_shade,
    shade_max               = maximum_shade,
    can_solar_orient        = true,
    can_press_to_ground     = true,
    can_change_absorptivity = true,
    alpha_min               = 0.7,
    alpha_max               = 0.9,
    alpha_step              = 0.01,
    can_pant                = false,
    pant_max                = 1.0,
    heat_exchange = example_ectotherm_heat_exchange_traits(;
        conduction_pars_external = example_ectotherm_conduction_pars_external(
            conduction_fraction = 0.1),
        evaporation_pars = example_ectotherm_evaporation_pars(
            eye_fraction = 0.0003, skin_wetness = 0.001),
        radiation_pars = example_ectotherm_radiation_pars(
            α_body_dorsal     = 0.85,
            α_body_ventral    = 0.85,
            solar_orientation = Intermediate(),
            ϵ_body_dorsal     = 0.95,
            ϵ_body_ventral    = 0.95),
        respiration_pars = example_ectotherm_respiration_pars(mouth_fraction = 0.0),
    ),
)

organism = Organism(body, organism_traits)
limits   = thermoregulation(organism)
env_pars = example_environment_pars(;
    elevation,
    α_ground = solar_terrain.albedo,
)

# ── Step 5: Thermoregulation loop ─────────────────────────────────────────
println("Running thermoregulation loop...")
results         = NamedTuple[]
prev_depth_node = limits.depth.reference
activity_today  = false

for step in 1:nsteps
    if (step - 1) % 24 == 0
        activity_today = false
    end
    out = thermoregulate(
        organism, available_environments, limits, env_pars, step, prev_depth_node;
        activity_today,
    )
    prev_depth_node = out.depth_node
    activity_today  = activity_today || out.state isa Active || out.state isa Basking
    push!(results, out)
end

# ── Extract outputs ───────────────────────────────────────────────────────
T_body   = [r.T_core     for r in results]
state    = [r.state      for r in results]
height   = [r.height     for r in results]
depth_nd = [r.depth_node for r in results]
shade    = [r.shade      for r in results]
T_air    = [low_shade_result.profile[i].air_temperature[1] for i in 1:nsteps]

T_body_C = ustrip.(u"°C", T_body)
act      = [s isa Active ? 2 : s isa Basking ? 1 : 0 for s in state]

_ground_ht = ustrip(u"cm", heights[1])
pos_cm = [depth_nd[i] > 1 ?
    -ustrip(u"cm", depths[depth_nd[i]]) :
    (h = ustrip(u"cm", height[i]); h > _ground_ht ? h : 0.0)
    for i in 1:nsteps]

month_ranges = [(m-1)*24+1 : m*24 for m in 1:ndays]
month_Tb     = [T_body_C[r]              for r in month_ranges]
month_act    = [act[r]                   for r in month_ranges]
month_Ta     = [ustrip.(u"°C", T_air[r]) for r in month_ranges]
month_pos    = [pos_cm[r]               for r in month_ranges]
month_shade  = [shade[r]                for r in month_ranges]

T_active_min_C = ustrip(u"°C", limits.T_active_min)
T_active_max_C = ustrip(u"°C", limits.T_active_max)

println("\n── Annual activity summary ──")
println("  Resting=$(sum(act.==0)), Basking=$(sum(act.==1)), Active=$(sum(act.==2))")

# ── Fig. 1 – Body temperature by month (4×3 grid) ────────────────────────
panels_Tb = map(1:ndays) do m
    p = plot(hours, month_Tb[m];
        lw = 2, color = :red, label = "",
        title = months[m], ylabel = "°C", ylim = (0, 50), titlefontsize = 9)
    plot!(p, hours, month_Ta[m]; lw = 1, color = :steelblue, linestyle = :dash, label = "")
    hline!(p, [T_active_min_C, T_active_max_C];
        color = :orange, linestyle = :dash, lw = 1, label = "")
    p
end

display(plot(panels_Tb...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Body temperature — generic lizard, Madison WI, 2000\n" *
                 "(red = Tb, blue dashed = T_air, orange dashed = active range)"))

# ── Fig. 2 – Annual heatmaps (body temp and activity) ─────────────────────
tb_matrix  = zeros(Float64, 24, ndays)
act_matrix = zeros(Int,     24, ndays)
for m in 1:ndays
    tb_matrix[:,  m] = month_Tb[m]
    act_matrix[:, m] = month_act[m]
end

p1 = heatmap(months, hours, tb_matrix;
    color = cgrad(:RdYlBu, rev = true), clims = (0, 50),
    colorbar_title = "°C",
    title = "Body temperature (°C)", ylabel = "hour")

p2 = heatmap(months, hours, act_matrix;
    color = cgrad([:steelblue, :orange, :firebrick], [0, 0.5, 1]),
    clims = (0, 2), colorbar_title = "0=rest 1=bask 2=active",
    title = "Activity state", ylabel = "hour")

display(plot(p1, p2; layout = (2, 1), size = (900, 600), left_margin = 6Plots.mm))

# ── Fig. 3 – Monthly activity budget ──────────────────────────────────────
n_rest   = [sum(month_act[m] .== 0) for m in 1:ndays]
n_bask   = [sum(month_act[m] .== 1) for m in 1:ndays]
n_active = [sum(month_act[m] .== 2) for m in 1:ndays]

display(groupedbar(hcat(n_rest, n_bask, n_active);
    bar_position = :stack,
    xticks       = (1:ndays, months),
    label        = ["rest" "bask" "active"],
    color        = [:steelblue :orange :firebrick],
    title        = "Monthly activity budget — generic lizard, Madison WI, 2000",
    ylabel       = "hours per day",
    ylim         = (0, 24),
    size         = (800, 400),
    left_margin  = 5Plots.mm,
))

# ── Fig. 4 – Shade selection by month (4×3 grid) ──────────────────────────
panels_shade = map(1:ndays) do m
    plot(hours, month_shade[m] .* 100;
        lw = 2, color = :darkgreen, label = "",
        title = months[m], ylabel = "%", ylim = (-5, 105), titlefontsize = 9)
end

display(plot(panels_shade...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Shade selection — generic lizard, Madison WI, 2000"))

# ── Fig. 5 – Position by month (height above / depth below ground) ─────────
depth_cm_max  = ustrip(u"cm", depths[end])
height_cm_max = ustrip(u"cm", heights[end])

panels_pos = map(1:ndays) do m
    p = plot(hours, month_pos[m];
        lw = 2, color = :sienna, label = "",
        title = months[m], ylabel = "cm",
        ylim = (-depth_cm_max, height_cm_max), titlefontsize = 9)
    hline!(p, [0.0]; color = :black, lw = 1, linestyle = :dash, label = "")
    hspan!(p, [-depth_cm_max, 0]; fillalpha = 0.06, color = :brown, label = "")
    hspan!(p, [0, height_cm_max]; fillalpha = 0.06, color = :skyblue, label = "")
    p
end

display(plot(panels_pos...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Position — generic lizard, Madison WI, 2000\n" *
                 "(above ground +cm, blue; underground −cm, brown)"))

# ── Fig. 6 – Annual heatmaps (shade and position) ─────────────────────────
# Colour map: stops are anchored to the actual discrete depth/height nodes so
# that underground always maps to brown, surface (0 cm) to limegreen, and
# above-ground to sky-blue — regardless of the depth/height vectors chosen.
pos_clims      = (-depth_cm_max, height_cm_max)
total_range    = depth_cm_max + height_cm_max
norm(v)        = (depth_cm_max + v) / total_range   # maps cm value → [0,1]
pos_shallowest = norm(-ustrip(u"cm", depths[2]))     # shallowest underground node
pos_surface    = norm(0.0)                           # surface
pos_lowest_ht  = norm(ustrip(u"cm", heights[2]))     # lowest above-ground node
pos_cmap = cgrad(
    [:saddlebrown, :saddlebrown, :limegreen, :limegreen, :skyblue, :steelblue],
    [0.0,
     (pos_shallowest + pos_surface) / 2,   # midpoint: shallowest underground → surface
     pos_surface - 0.001,
     pos_surface + 0.001,
     (pos_surface + pos_lowest_ht) / 2,    # midpoint: surface → lowest above-ground node
     1.0],
)

shade_matrix = zeros(Float64, 24, ndays)
pos_matrix   = zeros(Float64, 24, ndays)
for m in 1:ndays
    shade_matrix[:, m] = month_shade[m] .* 100
    pos_matrix[:,   m] = month_pos[m]
end

p_shade = heatmap(months, hours, shade_matrix;
    color = :Greens, clims = (0, 100),
    colorbar_title = "%", title = "Shade selection (%)", ylabel = "hour")

p_pos = heatmap(months, hours, pos_matrix;
    color = pos_cmap, clims = pos_clims,
    colorbar_title = "cm (+ above, − below)",
    title = "Position (cm above/below ground)", ylabel = "hour")

display(plot(p_shade, p_pos; layout = (2, 1), size = (900, 600), left_margin = 6Plots.mm))
