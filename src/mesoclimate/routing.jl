# Lateral surface-water routing along the D8 flow graph (flow_graph.jl). Cells solve in
# dependency order, so a cell's inflow is its upstream cells' already-computed
# `runoff_generated`, fed to the inner solve as a forcing. An atomic dataflow scheduler
# walks the dependency DAG; determinism holds across threads (fixed-order inflow sums).

abstract type RoutingModel end

"""
    SurfaceRunoffRouting(; method = D8(), sinks = TerminalSinks())

Route surface water (`runoff_generated`, the pool overflow above `max_surface_pool`)
downslope along the DEM flow graph, feeding it into each downstream column's water
balance as lateral inflow.

- `method` — only `D8()` (single steepest receiver) is supported.
- `sinks` — how closed depressions are treated: `TerminalSinks()` (default,
  endorheic accumulation) or `SpillOver()` (priority-flood, drains through).
"""
@kwdef struct SurfaceRunoffRouting{M,SH<:SinkHandling} <: RoutingModel
    method::M = D8()
    sinks::SH = TerminalSinks()
end

# Built once at `init` (grid mode only): the flow graph plus the reusable per-solve
# workspace — `store` (nsteps × ncells outflow, kg/m^2 stripped), one inflow `buffer`
# per worker, and the per-cell `remaining` upstream counters. `solve!` resets these
# rather than reallocating.
struct RoutingState{FG,M<:RoutingModel,ST,BUF}
    graph::FG
    model::M
    store::ST
    buffers::BUF
    remaining::Vector{Threads.Atomic{Int}}
    nsteps::Int
    n_active::Int
end

function build_routing_state(model::SurfaceRunoffRouting, elevation, mask, nsteps, nworkers)
    graph = build_flow_graph(elevation, mask; method = model.method, sinks = model.sinks)
    store = zeros(Float64, nsteps, graph.ncells)
    buffers = [zeros(typeof(0.0u"kg/m^2"), nsteps) for _ in 1:nworkers]
    remaining = [Threads.Atomic{Int}(0) for _ in 1:graph.ncells]
    return RoutingState(graph, model, store, buffers, remaining, nsteps, count(graph.active))
end

# Per-cell `max_surface_pool` by sink model; TerminalSinks gives sink cells its large
# pool so inflow accumulates there. New SinkHandling types add a method, not a branch.
_cell_pool(::SpillOver, _g, _id, normal) = normal
_cell_pool(s::TerminalSinks, g, id, normal) = g.is_sink[id] ? s.pool : normal

canonical_unit(::Val{:runoff_generated}) = u"kg/m^2"

# Sum upstream cells' stored outflow into the worker's buffer — this cell's inflow.
@inline function _gather_inflow!(buf, store, g::FlowGraph, id::Int)
    fill!(buf, 0.0u"kg/m^2")
    for k in _upstream_range(g, id)
        up = g.upstream_ids[k]
        for t in axes(store, 1)
            buf[t] += store[t, up] * u"kg/m^2"
        end
    end
    return buf
end

# Shared state for one routed solve; `store` is nsteps × ncells outflow (kg/m^2).
# One routed solve: the init-allocated `state` (graph + reusable workspace) plus the
# fields specific to this call — the output rasters, the ready queue, and progress.
struct RoutedSolve{C,S<:RoutingState,O,SO,RO,P}
    cache::C
    state::S
    output::O
    solar_output::SO
    runoff_out::RO
    normal_pool::P
    ready::Channel{Int}
    solved::Threads.Atomic{Int}
    show_progress::Bool
    t_start::Float64
    last_report_ms::Threads.Atomic{Int64}
end

@inline _has_solar(rs::RoutedSolve) = rs.solar_output !== nothing

# Clear the reused workspace so a fresh solve isn't contaminated by the previous one.
function _reset_workspace!(state::RoutingState)
    fill!(state.store, 0.0)
    for id in eachindex(state.remaining)
        state.remaining[id][] = state.graph.indegree[id]
    end
    return state
end

# Reset the init-allocated workspace for one solve and bind it to this call's outputs.
function _init_routed_solve(output, solar_output, cache, proto)
    state = cache.routing
    g = state.graph
    normal_pool = cache.problem.model.micro_model.config.max_surface_pool

    _reset_workspace!(state)
    ready = Channel{Int}(g.ncells)
    for id in 1:g.ncells
        g.active[id] && g.indegree[id] == 0 && put!(ready, id)   # in-degree-0 sources
    end
    ti = dims(first(values(output)), Ti)
    runoff_out = _allocate_layer_array(typeof(1.0u"kg/m^2"),
        (dims(cache.terrain.elevation)..., ti), cache.mask)

    t_start = time()
    return RoutedSolve(cache, state, output, solar_output, runoff_out, normal_pool,
        ready, Threads.Atomic{Int}(0), state.n_active > 10, t_start,
        Threads.Atomic{Int64}(round(Int64, (t_start - 10.1) * 1000)))
end

# Routing active: solve cells in flow-graph dependency order, routing runoff downslope.
function _solve_remaining!(::RoutingState, output, solar_output, cache, proto, _first_I)
    put!(cache.cache_pool, proto)          # proto is re-solved as an ordinary cell
    rs = _init_routed_solve(output, solar_output, cache, proto)

    @sync for w in 1:length(rs.state.buffers)
        Threads.@spawn begin
            c = take!(cache.cache_pool)
            try
                _run_worker!(rs, c, rs.state.buffers[w])
            finally
                put!(cache.cache_pool, c)
            end
        end
    end

    base = merge(NamedTuple(rs.output), (; runoff_generated = rs.runoff_out))
    _has_solar(rs) && (base = merge(base, NamedTuple(rs.solar_output)))
    return RasterStack(base)
end

# Walk each ready cell's flow path downstream, adopting cells we make ready.
function _run_worker!(rs::RoutedSolve, c, buf)
    for id in rs.ready
        while id != 0
            _solve_cell!(rs, c, buf, id)
            next_id = _release_downstream!(rs, id)
            _mark_solved!(rs)
            id = next_id
        end
    end
    return nothing
end

function _solve_cell!(rs::RoutedSolve, c, buf, id::Int)
    g = rs.state.graph
    I = g.dimindices[id]
    _gather_inflow!(buf, rs.state.store, g, id)
    pool = _cell_pool(rs.state.model.sinks, g, id, rs.normal_pool)
    reinit!(c.micro, rs.cache.init_inputs.build_inputs(c.scratch, I;
        lateral_inflow = buf, max_surface_pool = pool))
    ok = true
    try
        solve!(c.micro)             # one failed cell must not abort the run
    catch
        ok = false
    end
    if ok
        out = c.micro.output
        _write_output!(rs.output, out, rs.cache.problem.model.output_layers, I)
        _write_slice!(rs.runoff_out, out.runoff_generated, I)
        rs.state.store[:, id] .= ustrip.(u"kg/m^2", out.runoff_generated)
        _write_solar!(rs, c, I)
    else
        rs.state.store[:, id] .= 0.0      # failed solve sheds no water
    end
    return nothing
end

function _write_solar!(rs::RoutedSolve, c, I)
    _has_solar(rs) || return nothing
    _compute_solar_for_pixel!(c.scratch, rs.cache.terrain, rs.cache.albedo_grid, I)
    _write_solar_output!(rs.solar_output, c.scratch.solar.out,
        rs.cache.init_inputs.solar_pairs,
        rs.cache.cloud_constants.solar_model.wavelengths, I)
    return nothing
end

@inline function _release_downstream!(rs::RoutedSolve, id::Int)
    g = rs.state.graph
    r = g.receiver[id]
    (r != 0 && g.active[r] && Threads.atomic_sub!(rs.state.remaining[r], 1) == 1) ? r : 0
end

@inline function _mark_solved!(rs::RoutedSolve)
    n = Threads.atomic_add!(rs.solved, 1) + 1
    n == rs.state.n_active && close(rs.ready)
    rs.show_progress && _maybe_report_routing(n, rs.state.n_active, rs.t_start, rs.last_report_ms)
    return nothing
end

# Throttled (10 s) progress line, matching the independent loop's cadence.
function _maybe_report_routing(n, n_active, t_start, last_report_ms)
    t_now = time()
    t_now_ms = round(Int64, t_now * 1000)
    prev_ms = last_report_ms[]
    if t_now_ms - prev_ms >= 10_000 &&
            Threads.atomic_cas!(last_report_ms, prev_ms, t_now_ms) == prev_ms
        elapsed = t_now - t_start
        pct = round(Int, 100 * n / n_active)
        eta_s = round(Int, elapsed / n * (n_active - n))
        @info "Routed solve: $n / $n_active cells ($pct%) — ETA $(eta_s)s"
    end
    return nothing
end
