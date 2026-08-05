using Test
using MicroclimateMapper
using MicroclimateMapper: build_flow_graph, FlowGraph, _upstream_range
using Rasters, Rasters.Lookups
using Rasters: X, Y
using Unitful

# Synthetic elevation + all-active mask from a value matrix (metres). Ascending
# X/Y so the raster's own orientation matches the canonical matrix; the
# descending-axis path is exercised separately below.
function _elev_mask(values::AbstractMatrix; xstep = 100.0, ystep = 100.0)
    nx, ny = size(values)
    xs = X(range(0.0; step = xstep, length = nx); sampling = Intervals(Center()))
    ys = Y(range(0.0; step = ystep, length = ny); sampling = Intervals(Center()))
    elevation = Raster(Float64.(values), (xs, ys)) .* u"m"
    mask = Raster(trues(nx, ny), (xs, ys))
    return elevation, mask
end

# Follow `receiver` from `id`; return the outlet it reaches (0), or -1 on a cycle.
function _trace_to_outlet(g::FlowGraph, id::Int)
    steps = 0
    while id != 0
        id = g.receiver[id]
        steps += 1
        steps > g.ncells && return -1   # cycle guard
    end
    return 0
end

# Position of each cell id in the topological order, for the precedence invariant.
_order_pos(g::FlowGraph) = (p = zeros(Int, g.ncells); for (k, id) in enumerate(g.order); p[id] = k; end; p)

# Structural invariants that must hold for any valid drainage forest.
function _check_structure(g::FlowGraph)
    # Topological order covers every active cell exactly (no cycles).
    @test length(g.order) == count(g.active)
    @test allunique(g.order)
    @test all(id -> g.active[id], g.order)

    # Receivers are 0 or active; every active cell reaches an outlet.
    for id in 1:g.ncells
        g.active[id] || continue
        r = g.receiver[id]
        @test r == 0 || g.active[r]
        @test _trace_to_outlet(g, id) == 0
    end

    # CSR consistency: offsets monotone, upstream_ids partitions non-outlet cells.
    @test g.offsets[1] == 1
    @test g.offsets[end] - 1 == length(g.upstream_ids)
    @test length(g.upstream_ids) == sum(g.indegree)
    for id in 1:g.ncells
        @test g.indegree[id] == length(_upstream_range(g, id))
        for k in _upstream_range(g, id)
            @test g.receiver[g.upstream_ids[k]] == id      # every contributor drains here
        end
    end
    # Each active non-outlet cell appears exactly once as someone's upstream.
    contributors = sort(g.upstream_ids)
    expected = sort([id for id in 1:g.ncells if g.active[id] && g.receiver[id] != 0])
    @test contributors == expected

    # Precedence: every contributor precedes its receiver in the topological order.
    pos = _order_pos(g)
    for id in 1:g.ncells
        g.active[id] || continue
        r = g.receiver[id]
        r != 0 && @test pos[id] < pos[r]
    end

    # Sinks are terminal and active.
    for id in 1:g.ncells
        g.is_sink[id] || continue
        @test g.active[id]
        @test g.receiver[id] == 0
    end
end

@testset "tilted plane: strictly downhill, mass-conserving" begin
    # Elevation increases with X → interior water drains toward lower X. No pits,
    # so every receiver must be strictly non-uphill (validates index mapping /
    # receiver decoding — an axis flip would send flow uphill and fail here).
    nx, ny = 8, 6
    values = [10.0 * (i - 1) for i in 1:nx, _ in 1:ny]
    elevation, mask = _elev_mask(values)
    g = build_flow_graph(elevation, mask)

    _check_structure(g)

    # Downhill everywhere (no depressions on a monotonic plane).
    for id in 1:g.ncells
        g.active[id] || continue
        r = g.receiver[id]
        r == 0 && continue
        e_id = ustrip(u"m", elevation[g.dimindices[id]...])
        e_r  = ustrip(u"m", elevation[g.dimindices[r]...])
        @test e_r <= e_id + 1e-9
    end

    # Route 1 unit generated per cell downstream in topological order; all mass
    # must arrive at the outlets (open-boundary conservation).
    acc = zeros(Float64, g.ncells)
    for id in 1:g.ncells; g.active[id] && (acc[id] = 1.0); end
    for id in g.order
        r = g.receiver[id]
        r != 0 && (acc[r] += acc[id])
    end
    total_out = sum(acc[id] for id in 1:g.ncells if g.active[id] && g.receiver[id] == 0)
    @test total_out ≈ count(g.active)
end

@testset "single bowl: interior minimum is an endorheic sink" begin
    nx, ny = 9, 9
    cx, cy = (nx + 1) / 2, (ny + 1) / 2
    values = [ (i - cx)^2 + (j - cy)^2 for i in 1:nx, j in 1:ny ]  # paraboloid bowl

    # TerminalSinks (default) — steepest descent drains inward to the central
    # minimum, which becomes the single terminal sink.
    elevation, mask = _elev_mask(values)
    g = build_flow_graph(elevation, mask)
    _check_structure(g)
    @test count(g.is_sink) == 1
    sink = only(findall(g.is_sink))
    @test g.receiver[sink] == 0
    @test g.dimindices[sink] == (X(5), Y(5))       # geometric centre of the 9×9 bowl
    @test g.indegree[sink] > 0                      # water drains into it

    # SpillOver — priority-flood spills the pit over its rim; no sink.
    g2 = build_flow_graph(elevation, mask; sinks = SpillOver())
    _check_structure(g2)
    @test !any(g2.is_sink)
end

@testset "twin-basin ridge: independent drainage, invariants hold" begin
    # A north-south ridge down the middle splits the grid into two basins draining
    # to opposite edges. We assert the general invariants (acyclic, reachable,
    # CSR/precedence consistent) rather than exact receivers.
    nx, ny = 10, 6
    mid = (nx + 1) / 2
    values = [ -abs(i - mid) * 10.0 for i in 1:nx, _ in 1:ny ]  # peak at ridge, low at E/W edges
    elevation, mask = _elev_mask(values)
    g = build_flow_graph(elevation, mask)
    _check_structure(g)
end

@testset "descending Y axis maps back correctly" begin
    # Same monotonic plane but with a descending Y lookup: the canonical matrix
    # flips Y internally, and dimindices must map back so downhill still holds.
    nx, ny = 6, 6
    values = [10.0 * (i - 1) for i in 1:nx, _ in 1:ny]
    xs = X(range(0.0; step = 100.0, length = nx); sampling = Intervals(Center()))
    ys = Y(range(500.0; step = -100.0, length = ny); sampling = Intervals(Center()))  # descending
    elevation = Raster(Float64.(values), (xs, ys)) .* u"m"
    mask = Raster(trues(nx, ny), (xs, ys))
    g = build_flow_graph(elevation, mask)

    _check_structure(g)
    for id in 1:g.ncells
        g.active[id] || continue
        r = g.receiver[id]
        r == 0 && continue
        e_id = ustrip(u"m", elevation[g.dimindices[id]...])
        e_r  = ustrip(u"m", elevation[g.dimindices[r]...])
        @test e_r <= e_id + 1e-9
    end
end

@testset "masked cells are excluded from the graph" begin
    nx, ny = 6, 6
    values = [10.0 * (i - 1) for i in 1:nx, _ in 1:ny]
    elevation, mask = _elev_mask(values)
    mask[X(3), Y(3)] = false   # punch a hole
    g = build_flow_graph(elevation, mask)

    hole = findfirst(!, g.active)
    @test hole !== nothing
    @test g.dimindices[hole] == (X(3), Y(3))
    @test !g.active[hole]
    @test g.receiver[hole] == 0
    @test g.indegree[hole] == 0                       # nothing drains into a masked cell
    @test !any(==(hole), g.upstream_ids)              # nor is it anyone's contributor
    _check_structure(g)
end
