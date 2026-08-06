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

struct RoutingState{FG,M<:RoutingModel}
    graph::FG
    model::M
end

# Grid-mode only; points mode has no neighbour topology and rejects routing before this.
build_routing_state(model::SurfaceRunoffRouting, elevation, mask) =
    RoutingState(build_flow_graph(elevation, mask; method = model.method, sinks = model.sinks), model)

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
struct RoutedSolve{C,G<:FlowGraph,O,S,R,P,SK<:SinkHandling,B}
    cache::C
    g::G
    output::O
    solar_output::S
    store::Matrix{Float64}
    runoff_out::R
    normal_pool::P
    sinks::SK
    buffers::B                       # one reusable inflow buffer per worker
    nsteps::Int
    n_active::Int
    remaining::Vector{Threads.Atomic{Int}}
    ready::Channel{Int}
    solved::Threads.Atomic{Int}
    show_progress::Bool
    t_start::Float64
    last_report_ms::Threads.Atomic{Int64}
end

@inline _has_solar(rs::RoutedSolve) = rs.solar_output !== nothing

function _init_routed_solve(output, solar_output, cache, proto)
    g = cache.routing.graph
    normal_pool = cache.problem.model.micro_model.config.max_surface_pool

    nsteps = length(proto.micro.output.runoff_generated)
    store = zeros(Float64, nsteps, g.ncells)
    ti = dims(first(values(output)), Ti)
    runoff_out = _allocate_layer_array(typeof(1.0u"kg/m^2"),
        (dims(cache.terrain.elevation)..., ti), cache.mask)
    buffers = [zeros(typeof(0.0u"kg/m^2"), nsteps) for _ in 1:cache.cache_pool.sz_max]

    remaining = [Threads.Atomic{Int}(g.indegree[id]) for id in 1:g.ncells]
    ready = Channel{Int}(g.ncells)
    for id in 1:g.ncells
        g.active[id] && g.indegree[id] == 0 && put!(ready, id)   # in-degree-0 sources
    end

    n_active = count(g.active)
    t_start = time()
    return RoutedSolve(cache, g, output, solar_output, store, runoff_out,
        normal_pool, cache.routing.model.sinks, buffers, nsteps, n_active,
        remaining, ready, Threads.Atomic{Int}(0), n_active > 10, t_start,
        Threads.Atomic{Int64}(round(Int64, (t_start - 10.1) * 1000)))
end

function _solve_routed!(output, solar_output, cache, proto)
    put!(cache.cache_pool, proto)          # proto is re-solved as an ordinary cell
    rs = _init_routed_solve(output, solar_output, cache, proto)

    @sync for w in 1:length(rs.buffers)
        Threads.@spawn begin
            c = take!(cache.cache_pool)
            try
                _run_worker!(rs, c, rs.buffers[w])
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
    g = rs.g
    I = g.dimindices[id]
    _gather_inflow!(buf, rs.store, g, id)
    pool = _cell_pool(rs.sinks, g, id, rs.normal_pool)
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
        rs.store[:, id] .= ustrip.(u"kg/m^2", out.runoff_generated)
        _write_solar!(rs, c, I)
    else
        rs.store[:, id] .= 0.0      # failed solve sheds no water
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
    r = rs.g.receiver[id]
    (r != 0 && rs.g.active[r] && Threads.atomic_sub!(rs.remaining[r], 1) == 1) ? r : 0
end

@inline function _mark_solved!(rs::RoutedSolve)
    n = Threads.atomic_add!(rs.solved, 1) + 1
    n == rs.n_active && close(rs.ready)
    rs.show_progress && _maybe_report_routing(n, rs.n_active, rs.t_start, rs.last_report_ms)
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
