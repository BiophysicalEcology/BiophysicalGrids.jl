using Test
using MicroclimateMapper
using MicroclimateMapper: build_routing_state, RoutingState, build_flow_graph,
    _gather_inflow!, _cell_pool, _reset_workspace!, _upstream_range
using Rasters, Rasters.Lookups
using Rasters: X, Y
using Unitful

# Synthetic elevation + all-active mask from a value matrix (metres).
function _elev_mask(values::AbstractMatrix; xstep = 100.0, ystep = 100.0)
    nx, ny = size(values)
    xs = X(range(0.0; step = xstep, length = nx); sampling = Intervals(Center()))
    ys = Y(range(0.0; step = ystep, length = ny); sampling = Intervals(Center()))
    elevation = Raster(Float64.(values), (xs, ys)) .* u"m"
    mask = Raster(trues(nx, ny), (xs, ys))
    return elevation, mask
end

# 3×3 bowl: the centre is a local minimum, so every edge cell drains inward and the
# centre becomes an interior endorheic sink (indegree 8, receiver 0).
_bowl() = _elev_mask([2.0 2.0 2.0; 2.0 1.0 2.0; 2.0 2.0 2.0])

@testset "build_routing_state preallocates the reusable workspace at init" begin
    elevation, mask = _bowl()
    nsteps, nworkers = 4, 3
    rs = build_routing_state(SurfaceRunoffRouting(), elevation, mask, nsteps, nworkers)

    @test rs isa RoutingState
    @test size(rs.store) == (nsteps, rs.graph.ncells)
    @test all(iszero, rs.store)
    @test length(rs.buffers) == nworkers
    @test all(b -> length(b) == nsteps, rs.buffers)
    @test length(rs.remaining) == rs.graph.ncells
    @test rs.nsteps == nsteps
    @test rs.n_active == count(rs.graph.active)
end

@testset "_reset_workspace! clears state so a reused cache isn't contaminated" begin
    elevation, mask = _bowl()
    rs = build_routing_state(SurfaceRunoffRouting(), elevation, mask, 4, 2)

    # Dirty the workspace the way a completed solve leaves it, then reset.
    rs.store .= 99.0
    for a in rs.remaining
        a[] = -1
    end
    _reset_workspace!(rs)

    @test all(iszero, rs.store)
    @test all(id -> rs.remaining[id][] == rs.graph.indegree[id], eachindex(rs.remaining))
end

@testset "_gather_inflow! sums upstream cells' stored outflow" begin
    elevation, mask = _bowl()
    nsteps = 3
    rs = build_routing_state(SurfaceRunoffRouting(), elevation, mask, nsteps, 1)
    g = rs.graph

    sink = findfirst(g.is_sink)
    ups = [g.upstream_ids[k] for k in _upstream_range(g, sink)]
    @test length(ups) == 8                         # all eight neighbours drain in
    for up in ups
        rs.store[:, up] .= Float64(up)             # distinct outflow per contributor
    end

    buf = zeros(typeof(0.0u"kg/m^2"), nsteps)
    _gather_inflow!(buf, rs.store, g, sink)
    @test all(t -> buf[t] ≈ sum(Float64.(ups)) * u"kg/m^2", 1:nsteps)

    # A source cell (indegree 0) gathers nothing, regardless of the store.
    src = findfirst(id -> g.indegree[id] == 0 && g.active[id], 1:g.ncells)
    _gather_inflow!(buf, rs.store, g, src)
    @test all(iszero, buf)
end

@testset "_cell_pool dispatches on the sink model" begin
    normal = 1.0u"kg/m^2"
    elevation, mask = _bowl()

    # SpillOver has no sinks — always the normal pool.
    g_spill = build_flow_graph(elevation, mask; sinks = SpillOver())
    @test _cell_pool(SpillOver(), g_spill, 1, normal) == normal

    # TerminalSinks — the sink cell gets the model's (large) pool, others the normal one.
    g_term = build_flow_graph(elevation, mask; sinks = TerminalSinks())
    sink = findfirst(g_term.is_sink)
    nonsink = findfirst(id -> !g_term.is_sink[id] && g_term.active[id], 1:g_term.ncells)
    s = TerminalSinks(pool = 5.0u"kg/m^2")
    @test _cell_pool(s, g_term, sink, normal) == s.pool
    @test _cell_pool(s, g_term, nonsink, normal) == normal
end
