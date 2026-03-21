# Endotherm thermoregulation at Madison, WI using TerraClimate year 2000.
#
# Demonstrates the full pipeline:
#   1. get_weather / apply_climate_scenario  — TerraClimate forcing (BiophysicalGrids)
#   2. simulate_microclimate                 — hourly air temp, humidity, wind (BiophysicalGrids)
#   3. thermoregulate                        — endotherm physiology (BiophysicalBehaviour)
#
# The organism uses BiophysicalBehaviour's default NicheMapR-equivalent parameters
# (~65 kg generic mammal with fur insulation).
#
# Requires BiophysicalBehaviour.jl (not a core BiophysicalGrids dependency):
#   using Pkg; Pkg.add(url="https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl")

using BiophysicalGrids
using BiophysicalBehaviour
using HeatExchange
using BiophysicalGeometry
using Microclimate
using SolarRadiation
using FluidProperties
using ConstructionBase
using ModelParameters
using Unitful, UnitfulMoles
using Statistics
using Plots

# ── Location ──────────────────────────────────────────────────────────────
lon, lat  = -89.4557, 43.1379
elevation = 270.0u"m"

months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
hours  = collect(0.0:1.0:23.0)
ndays  = 12
nsteps = ndays * 24

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

# ── Step 3: Solve microclimate ────────────────────────────────────────────
# Heights [0.01, 2.0] m: index 1 = near-ground, index 2 = 2 m reference height.
println("Solving microclimate...")
micro_result = simulate_microclimate(
    solar_terrain, micro_terrain, soil_thermal_model, weather_scenario;
    depths   = [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]u"cm",
    heights  = [0.01, 2.0]u"m",
    runmoist = false,
    organic_soil_cap = true,
)

# ── Step 4: Set up endotherm ──────────────────────────────────────────────
# Default NicheMapR-equivalent parameters for a generic ~65 kg mammal with fur.
shape_pars       = example_shape_pars()
insulation_pars  = example_insulation_pars(; 
                    insulation_depth_dorsal = 2.0u"mm",
                    insulation_depth_ventral = 2.0u"mm",
                    )
radiation_pars   = example_radiation_pars()
metabolism_pars  = example_metabolism_pars()
evaporation_pars = example_evaporation_pars()
respiration_pars = example_respiration_pars()

conduction_pars_internal = example_conduction_pars_internal()
fat = Fat(conduction_pars_internal.fat_fraction, conduction_pars_internal.ρ_fat)

mean_insulation_depth = insulation_pars.dorsal.depth * (1 - radiation_pars.ventral_fraction) +
                        insulation_pars.ventral.depth * radiation_pars.ventral_fraction
mean_fibre_diameter   = insulation_pars.dorsal.diameter * (1 - radiation_pars.ventral_fraction) +
                        insulation_pars.ventral.diameter * radiation_pars.ventral_fraction
mean_fibre_density    = insulation_pars.dorsal.density * (1 - radiation_pars.ventral_fraction) +
                        insulation_pars.ventral.density * radiation_pars.ventral_fraction
fur      = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

# Plot the initial shape if interested
# using CairoMakie
# fig = plot_cross_sections(geometry)

physiology_traits = HeatExchangeTraits(
    shape_pars,
    insulation_pars,
    example_conduction_pars_external(),
    conduction_pars_internal,
    radiation_pars,
    ConvectionParameters(),
    evaporation_pars,
    example_hydraulic_pars(),
    respiration_pars,
    metabolism_pars,
    example_metabolic_rate_options(),
)

T_core_ref = metabolism_pars.T_core

thermoregulation_limits = ThermoregulationLimits(;
    control       = RuleBasedSequentialControl(; mode = 1, tolerance = 0.005, max_iterations = 200),
    Q_minimum_ref = metabolism_pars.Q_metabolism,
    insulation    = InsulationLimits(;
        dorsal  = SteppedParameter(; current   = insulation_pars.dorsal.depth,
                                     reference = insulation_pars.dorsal.depth,
                                     max       = insulation_pars.dorsal.depth, step = 0.0),
        ventral = SteppedParameter(; current   = insulation_pars.ventral.depth,
                                     reference = insulation_pars.ventral.depth,
                                     max       = insulation_pars.ventral.depth, step = 0.0),
    ),
    shape_b      = SteppedParameter(; current = shape_pars.b, max = 5.0, step = 0.1),
    k_flesh      = SteppedParameter(; current = conduction_pars_internal.k_flesh,
                                       max = 2.8u"W/m/K", step = 0.1u"W/m/K"),
    T_core       = SteppedParameter(; current = T_core_ref, reference = T_core_ref,
                                       max = T_core_ref + 5.0u"K", step = 0.1u"K"),
    panting      = PantingLimits(;
        pant       = SteppedParameter(; current = respiration_pars.pant, max = 15.0, step = 0.01),
        cost       = 0.0u"W",
        multiplier = 1.0,
        T_core_ref,
    ),
    skin_wetness = SteppedParameter(; current = evaporation_pars.skin_wetness,
                                       max = 0.5, step = 0.01),
)

behavioral_traits = BehavioralTraits(;
    thermoregulation = thermoregulation_limits,
    activity_period  = Diurnal(),
)
organism_traits = OrganismTraits(Endotherm(), physiology_traits, behavioral_traits)
organism        = Organism(geometry, organism_traits)

environment_pars = example_environment_pars(; elevation)

# ── Step 5: Thermoregulation loop ─────────────────────────────────────────
# Warm-start: carry T_skin and T_insulation forward between hours for faster convergence.
T_skin_prev       = T_core_ref - 3.0u"K"
T_insulation_prev = u"K"(10.0u"°C")
Q_gen             = 0.0u"W"

endo_results = Vector{NamedTuple}(undef, nsteps)
println("Running endotherm thermoregulation loop...")

for step in 1:nsteps
    profile = micro_result.profile[step]
    T_air   = profile.air_temperature[2]    # 2 m reference height
    rh      = profile.relative_humidity[2]
    wind    = profile.wind_speed[2]

    environment_vars = example_environment_vars(;
        T_air,
        rh,
        wind_speed       = wind,
        P_atmos          = atmospheric_pressure(elevation),
        global_radiation = 0.0u"W/m^2",   # endotherm assumed in shade/shelter
        zenith_angle     = 20.0u"°",
    )

    out = thermoregulate(
        organism,
        (; environment_pars, environment_vars),
        Q_gen,
        T_skin_prev,
        T_insulation_prev,
    )

    endo_results[step]  = out
    T_skin_prev         = out.thermoregulation.T_skin
    T_insulation_prev   = out.thermoregulation.T_insulation
end

# ── Extract outputs ───────────────────────────────────────────────────────
T_air_C      = [ustrip(u"°C",   micro_result.profile[i].air_temperature[2])  for i in 1:nsteps]
T_core_C     = [ustrip(u"°C",   endo_results[i].thermoregulation.T_core)     for i in 1:nsteps]
Q_gen_W      = [ustrip(u"W",    endo_results[i].energy_fluxes.Q_gen)         for i in 1:nsteps]
m_evap_gh    = [ustrip(u"g/hr", endo_results[i].mass_fluxes.m_evap)          for i in 1:nsteps]
m_resp_gh    = [ustrip(u"g/hr", endo_results[i].mass_fluxes.m_resp)          for i in 1:nsteps]
m_sweat_gh   = [ustrip(u"g/hr", endo_results[i].mass_fluxes.m_sweat)         for i in 1:nsteps]
shape_b      = [endo_results[i].thermoregulation.shape_b                      for i in 1:nsteps]
skin_wetness = [endo_results[i].thermoregulation.skin_wetness                 for i in 1:nsteps]
pant         = [endo_results[i].thermoregulation.pant                         for i in 1:nsteps]

month_ranges    = [(m-1)*24+1 : m*24 for m in 1:ndays]
month_Ta        = [T_air_C[r]       for r in month_ranges]
month_Tc        = [T_core_C[r]      for r in month_ranges]
month_Qg        = [Q_gen_W[r]       for r in month_ranges]
month_evap      = [m_evap_gh[r]     for r in month_ranges]
month_resp      = [m_resp_gh[r]     for r in month_ranges]
month_sweat     = [m_sweat_gh[r]    for r in month_ranges]
month_shape_b   = [shape_b[r]       for r in month_ranges]
month_wetness   = [skin_wetness[r]  for r in month_ranges]
month_pant      = [pant[r]          for r in month_ranges]

println("\n── Annual metabolic summary ──")
println("  Mean Q_gen: $(round(mean(Q_gen_W); digits=2)) W")
println("  Max  Q_gen: $(round(maximum(Q_gen_W); digits=2)) W")
println("  Min  Q_gen: $(round(minimum(Q_gen_W); digits=2)) W")
println("\n── Annual water loss summary ──")
println("  Mean total evap: $(round(mean(m_evap_gh); digits=3)) g/hr")
println("  Mean resp loss:  $(round(mean(m_resp_gh); digits=3)) g/hr")
println("  Mean cutaneous:  $(round(mean(m_sweat_gh); digits=3)) g/hr")

# ── Fig. 1 – Core temperature by month (4×3 grid) ────────────────────────
panels_Tc = map(1:ndays) do m
    p = plot(hours, month_Tc[m];
        lw = 2, color = :red, label = "",
        title = months[m], ylabel = "°C", titlefontsize = 9)
    plot!(p, hours, month_Ta[m];
        lw = 1, color = :steelblue, linestyle = :dash, label = "")
    p
end

display(plot(panels_Tc...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Core temperature — generic endotherm, Madison WI, 2000\n" *
                 "(red = T_core, blue dashed = T_air)"))

# ── Fig. 2 – Metabolic heat by month (4×3 grid) ───────────────────────────
panels_Qg = map(1:ndays) do m
    plot(hours, month_Qg[m];
        lw = 2, color = :orange, label = "",
        title = months[m], ylabel = "W", titlefontsize = 9)
end

display(plot(panels_Qg...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Metabolic heat (Q_gen) — generic endotherm, Madison WI, 2000"))

# ── Fig. 3 – Annual heatmaps (T_core and Q_gen) ───────────────────────────
tc_matrix = zeros(Float64, 24, ndays)
qg_matrix = zeros(Float64, 24, ndays)
for m in 1:ndays
    tc_matrix[:, m] = month_Tc[m]
    qg_matrix[:, m] = month_Qg[m]
end

p_tc = heatmap(months, hours, tc_matrix;
    color = cgrad(:RdYlBu, rev = true),
    colorbar_title = "°C",
    title = "Core temperature (°C)", ylabel = "hour")

p_qg = heatmap(months, hours, qg_matrix;
    color = :heat,
    colorbar_title = "W",
    title = "Metabolic heat Q_gen (W)", ylabel = "hour")

display(plot(p_tc, p_qg; layout = (2, 1), size = (900, 600), left_margin = 6Plots.mm))

# ── Fig. 4 – Water loss by month (4×3 grid) ───────────────────────────────
panels_wl = map(1:ndays) do m
    p = plot(hours, month_evap[m];
        lw = 2, color = :teal, label = "total",
        title = months[m], ylabel = "g/hr", titlefontsize = 9)
    plot!(p, hours, month_resp[m];
        lw = 1, color = :purple, linestyle = :dash, label = "resp")
    plot!(p, hours, month_sweat[m];
        lw = 1, color = :orange, linestyle = :dot, label = "cutaneous")
    p
end

display(plot(panels_wl...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Water loss — generic endotherm, Madison WI, 2000\n" *
                 "(teal = total, purple dashed = respiratory, orange dotted = cutaneous)"))

# ── Fig. 5 – Annual heatmap of total water loss ───────────────────────────
evap_matrix = zeros(Float64, 24, ndays)
for m in 1:ndays
    evap_matrix[:, m] = month_evap[m]
end

p_evap = heatmap(months, hours, evap_matrix;
    color = :Blues,
    colorbar_title = "g/hr",
    title = "Total evaporative water loss (g/hr)", ylabel = "hour")

display(plot(p_evap; size = (900, 350), left_margin = 6Plots.mm))

# ── Fig. 6 – Posture (shape_b) by month (4×3 grid) ────────────────────────
sb_ylim = (0.9, max(5.1, maximum(vcat(month_shape_b...)) + 0.2))
panels_sb = map(1:ndays) do m
    plot(hours, month_shape_b[m];
        lw = 2, color = :sienna, label = "",
        title = months[m], ylabel = "shape b",
        ylim = sb_ylim, titlefontsize = 9)
end

display(plot(panels_sb...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Posture (shape b) — generic endotherm, Madison WI, 2000"))

# ── Fig. 7 – Skin wetness by month (4×3 grid) ────────────────────────────
panels_sw = map(1:ndays) do m
    plot(hours, month_wetness[m];
        lw = 2, color = :dodgerblue, label = "",
        title = months[m], ylabel = "skin wetness",
        ylim = (0.0, max(0.05, maximum(vcat(month_wetness...)) * 1.1)),
        titlefontsize = 9)
end

display(plot(panels_sw...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Skin wetness — generic endotherm, Madison WI, 2000"))

# ── Fig. 8 – Pant rate by month (4×3 grid) ───────────────────────────────
pt_ylim = (0.9, max(1.1, maximum(vcat(month_pant...)) * 1.05))
panels_pt = map(1:ndays) do m
    plot(hours, month_pant[m];
        lw = 2, color = :crimson, label = "",
        title = months[m], ylabel = "pant rate",
        ylim = pt_ylim, titlefontsize = 9)
end

display(plot(panels_pt...; layout = (4, 3), size = (1200, 900),
    xlabel = "hour", left_margin = 4Plots.mm,
    plot_title = "Pant rate — generic endotherm, Madison WI, 2000"))

# ── Fig. 9 – Annual heatmaps (posture, skin wetness, panting) ─────────────
sb_matrix  = zeros(Float64, 24, ndays)
sw_matrix  = zeros(Float64, 24, ndays)
pt_matrix  = zeros(Float64, 24, ndays)
for m in 1:ndays
    sb_matrix[:, m] = month_shape_b[m]
    sw_matrix[:, m] = month_wetness[m]
    pt_matrix[:, m] = month_pant[m]
end

p_sb = heatmap(months, hours, sb_matrix;
    color = :YlOrBr, colorbar_title = "",
    title = "Posture (shape b)", ylabel = "hour")
p_sw = heatmap(months, hours, sw_matrix;
    color = :Blues, colorbar_title = "",
    title = "Skin wetness", ylabel = "hour")
p_pt = heatmap(months, hours, pt_matrix;
    color = :Reds, colorbar_title = "",
    title = "Pant rate", ylabel = "hour")

display(plot(p_sb, p_sw, p_pt; layout = (3, 1), size = (900, 800), left_margin = 6Plots.mm))
