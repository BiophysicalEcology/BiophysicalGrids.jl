# routing.jl — lateral surface-water routing along the D8 flow graph (flow_graph.jl).
# Cells solve in dependency order, so a cell's lateral inflow is the already-computed
# `runoff_generated` series of its upstream cells, fed to the inner solve as a forcing
# (`MicroInputs.lateral_inflow`); its overflow routes further downslope. An atomic
# dataflow scheduler exploits the width of the dependency DAG — workers seed from
# in-degree-0 sources and release downstream cells via atomic in-degree counters.
# Determinism holds across thread counts: each inflow is a fixed-order serial sum.

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

# Per-run workspace: flow graph + model. The per-cell time-series store is
# allocated inside the solve, so this is cheap to build at `init`.
struct RoutingState{FG,M<:RoutingModel}
    graph::FG
    model::M
end

# Build the flow graph at `init` time and pair it with the routing model. Grid-mode
# only — points mode has no neighbour topology and rejects routing before this.
build_routing_state(model::SurfaceRunoffRouting, elevation, mask) =
    RoutingState(build_flow_graph(elevation, mask; method = model.method, sinks = model.sinks), model)

# The `max_surface_pool` used for sink cells: TerminalSinks carries it; SpillOver
# has no sinks, so this is never applied (fall back to the model's normal pool).
_sink_pool(s::TerminalSinks, _normal) = s.pool
_sink_pool(::SpillOver, normal) = normal

# Canonical reporting unit for the routed-runoff output layer.
canonical_unit(::Val{:runoff_generated}) = u"kg/m^2"

# This cell's lateral inflow (kg/m^2 per output step): sum of upstream cells'
# stored outflow, into the worker's reusable buffer. No transmission loss — the
# column's own water balance removes what doesn't run on. Fixed order ⇒
# thread-count-independent.
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

# ---------------------------------------------------------------------------
# Atomic dataflow scheduler
# ---------------------------------------------------------------------------

# Routed variant of `_solve_remaining!`. Solves every active cell in dependency
# order across the worker pool, threading each cell's upstream outflow in as
# `lateral_inflow` and recording its own `runoff_generated` for downstream cells.
function _solve_routed!(output, solar_output, cache, proto, first_I)
    rs = cache.routing
    g = rs.graph
    model = cache.problem.model
    # Endorheic sink cells get an effectively unbounded pool so routed inflow
    # accumulates (removed only by evaporation/infiltration); all other cells
    # keep the model's normal max_surface_pool and shed overflow downslope.
    normal_pool = model.micro_model.config.max_surface_pool
    sink_pool = _sink_pool(rs.model.sinks, normal_pool)
    layers = model.output_layers
    build_inputs = cache.init_inputs.build_inputs
    has_solar = solar_output !== nothing
    solar_pairs = cache.init_inputs.solar_pairs
    wavelengths = has_solar ? cache.cloud_constants.solar_model.wavelengths : nothing

    nsteps = length(proto.micro.output.runoff_generated)
    put!(cache.cache_pool, proto)          # proto is re-solved as an ordinary cell

    ncells = g.ncells
    n_active = count(g.active)
    # Per-cell lateral runoff time series (kg/m^2, stripped), column = one cell.
    store = zeros(Float64, nsteps, ncells)
    # Runoff output raster (merged into the returned stack), shaped like a scalar layer.
    ti = dims(first(values(output)), Ti)
    runoff_out = _allocate_layer_array(typeof(1.0u"kg/m^2"),
        (dims(cache.terrain.elevation)..., ti), cache.mask)

    # Remaining upstream-dependency counters; ready = in-degree-0 sources.
    remaining = [Threads.Atomic{Int}(g.indegree[id]) for id in 1:ncells]
    ready = Channel{Int}(ncells)
    for id in 1:ncells
        g.active[id] && g.indegree[id] == 0 && put!(ready, id)
    end

    solved = Threads.Atomic{Int}(0)
    show_progress = n_active > 10
    t_start = time()
    last_report_ms = Threads.Atomic{Int64}(round(Int64, (t_start - 10.1) * 1000))

    nworkers = cache.cache_pool.sz_max
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            c = take!(cache.cache_pool)
            buf = c.scratch.lateral_inflow          # worker-owned, reused every cell
            try
                while true
                    id = try
                        take!(ready)
                    catch e
                        e isa InvalidStateException && break   # queue closed → done
                        rethrow()
                    end
                    # Walk down the flow path: adopt each downstream cell that this
                    # worker makes ready (cache-friendly; no queue round-trip).
                    while id != 0
                        I = g.dimindices[id]
                        _gather_inflow!(buf, store, g, id)
                        c.scratch.max_surface_pool[] = g.is_sink[id] ? sink_pool : normal_pool
                        reinit!(c.micro, build_inputs(c.scratch, I))
                        ok = true
                        try
                            solve!(c.micro)
                        catch
                            ok = false
                        end
                        if ok
                            _write_output!(output, c.micro.output, layers, I)
                            _write_slice!(runoff_out, c.micro.output.runoff_generated, I)
                            for t in 1:nsteps
                                store[t, id] = ustrip(u"kg/m^2", c.micro.output.runoff_generated[t])
                            end
                            if has_solar
                                _compute_solar_for_pixel!(c.scratch, cache.terrain, cache.albedo_grid, I)
                                _write_solar_output!(solar_output, c.scratch.solar.out,
                                    solar_pairs, wavelengths, I)
                            end
                        else
                            for t in 1:nsteps
                                store[t, id] = 0.0   # failed solve sheds no water
                            end
                        end

                        # Release the single downstream cell; adopt it if we made it ready.
                        next_id = 0
                        r = g.receiver[id]
                        if r != 0 && g.active[r] && Threads.atomic_sub!(remaining[r], 1) == 1
                            next_id = r
                        end

                        n = Threads.atomic_add!(solved, 1) + 1
                        n == n_active && close(ready)
                        if show_progress
                            _maybe_report_routing(n, n_active, t_start, last_report_ms)
                        end
                        id = next_id
                    end
                end
            finally
                put!(cache.cache_pool, c)
            end
        end
    end

    base = merge(NamedTuple(output), (; runoff_generated = runoff_out))
    has_solar && (base = merge(base, NamedTuple(solar_output)))
    return RasterStack(base)
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
