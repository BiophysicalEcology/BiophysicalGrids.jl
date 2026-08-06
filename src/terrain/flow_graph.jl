# flow_graph.jl — D8 drainage forest built once from the run-grid DEM. Each active
# cell drains to one downslope neighbour, so cells solve in dependency order and a
# downstream cell's inflow is its upstream cells' already-computed outflow series
# (a forcing, no global time-stepping). Stored as flat arrays for the routing.jl
# scheduler: `receiver` (forward edge) + a CSR pair `upstream_ids`/`offsets`.

# Directed D8 drainage forest over the run grid. Ids are linear indices over the
# `elevation` grid; `dimindices` maps each id back to its raster selector for
# `terrain`/`output`/`build_inputs`.
# Fields: `receiver` downstream id per cell (0 = outlet/inactive); `upstream_ids`
# + `offsets` a CSR store of contributors (cell `id`'s are
# `upstream_ids[offsets[id]:offsets[id+1]-1]`); `indegree` contributor count;
# `order` one topological order (sources first); `active`/`is_sink` bit masks —
# `is_sink` marks interior endorheic minima, distinct from domain-edge outlets.
struct FlowGraph{DI}
    receiver::Vector{Int}
    upstream_ids::Vector{Int}
    offsets::Vector{Int}
    indegree::Vector{Int}
    order::Vector{Int}
    dimindices::DI
    active::BitVector
    is_sink::BitVector
    ncells::Int
end

# 8-neighbour offsets and their in-cell distances (orthogonal 1, diagonal √2),
# for steepest-descent flow directions.
const _NB8 = (
    CartesianIndex(-1, 0), CartesianIndex(1, 0), CartesianIndex(0, -1), CartesianIndex(0, 1),
    CartesianIndex(-1, -1), CartesianIndex(-1, 1), CartesianIndex(1, -1), CartesianIndex(1, 1),
)
const _NB8_DIST = (1.0, 1.0, 1.0, 1.0, sqrt(2.0), sqrt(2.0), sqrt(2.0), sqrt(2.0))

# Index range of cell `id`'s upstream contributors in the CSR `upstream_ids` store.
@inline _upstream_range(g::FlowGraph, id::Int) = g.offsets[id]:(g.offsets[id+1] - 1)

# How closed topographic depressions are treated when the graph is built —
# whether water at a local minimum retains (`TerminalSinks`) or spills (`SpillOver`).
abstract type SinkHandling end

"""
    TerminalSinks(; pool = 1.0e12u"kg/m^2")

Closed depressions are terminal endorheic sinks: routed water accumulates in the
local minimum and leaves only vertically via the column's own water balance, never
laterally. Steepest-descent flow, so undrained flats are terminal too. Default.

- `pool` — `max_surface_pool` for sink cells; large enough that inflow accumulates
  and is drawn down only by evaporation and infiltration.
"""
@kwdef struct TerminalSinks <: SinkHandling
    pool::typeof(1.0u"kg/m^2") = 1.0e12u"kg/m^2"
end

"""
    SpillOver()

Closed depressions are assumed already full and overflowing at their sill, passing
all inflow downslope with no net storage — an exorheic network draining to the
domain edge everywhere. Flow directions come from a priority-flood over every pit.
"""
struct SpillOver <: SinkHandling end

# Build the D8 drainage forest from the run-grid `elevation` (Unitful m) and the
# active `mask` (2-D Bool on the same grid); masked/non-finite cells never become a
# receiver. `sinks` picks depression handling (TerminalSinks endorheic default, or
# SpillOver priority-flood) — both strictly downhill and acyclic, so the scheduler
# and order are identical. Only single-receiver `D8()` is supported.
function build_flow_graph(elevation::Raster, mask; method = D8(), sinks::SinkHandling = TerminalSinks())
    method isa D8 || error(
        "build_flow_graph: only D8() routing is supported (single-receiver flow " *
        "with reconstructable weights); got $(method).")

    # elevation (Float64, stripped) + active mask, addressed by `DimIndices` so raster
    # storage orientation is irrelevant — D8 drainage is orientation-invariant.
    dimidx = DimIndices(elevation)
    dem = zeros(dims(dimidx))
    active2d = falses(dims(dimidx))
    for D in dimidx
        v = ustrip(u"m", elevation[D])
        on = isfinite(v) && Bool(mask[D])
        active2d[D] = on
        dem[D] = on ? v : 0.0   # value under masked cells is never read by the flood
    end

    R = CartesianIndices(dem)
    L = LinearIndices(dem)
    ncells = length(dem)
    active = BitVector(vec(active2d))

    receiver = _receivers(sinks, dem, active2d, R, L)

    # Terminal cells ringed entirely by active in-grid cells are endorheic sinks;
    # those touching the domain edge or a masked region are outlets.
    is_sink = falses(ncells)
    for c in R
        active2d[c] || continue
        id = L[c]
        receiver[id] == 0 || continue
        interior = true
        for off in _NB8
            nc = c + off
            if !(nc in R) || !active2d[nc]
                interior = false
                break
            end
        end
        is_sink[id] = interior
    end

    # CSR inverse adjacency (upstream contributors), built by counting sort.
    indegree = zeros(Int, ncells)
    for id in 1:ncells
        r = receiver[id]
        r == 0 && continue
        indegree[r] += 1
    end
    offsets = Vector{Int}(undef, ncells + 1)
    offsets[1] = 1
    for id in 1:ncells
        offsets[id + 1] = offsets[id] + indegree[id]
    end
    upstream_ids = Vector{Int}(undef, offsets[end] - 1)
    cursor = copy(offsets)
    for id in 1:ncells
        r = receiver[id]
        r == 0 && continue
        upstream_ids[cursor[r]] = id
        cursor[r] += 1
    end

    order = _topological_order(receiver, indegree, active, ncells)
    dimindices = collect(vec(dimidx))   # linear cell id → raster dim-index

    return FlowGraph(receiver, upstream_ids, offsets, indegree, order,
                     dimindices, active, is_sink, ncells)
end

# Flow directions per sink-handling model.
_receivers(::TerminalSinks, dem, active2d, R, L) = _steepest_descent_receivers(dem, active2d, R, L)
_receivers(::SpillOver, dem, active2d, R, L) = _priorityflood_receivers(dem, active2d, R, L)

# Steepest-descent D8: each active cell drains to its strictly-lowest in-grid
# active neighbour (drop / distance, orthogonal 1, diagonal √2). A cell with no
# lower neighbour has `receiver = 0` — terminal (a pit/flat that holds water).
function _steepest_descent_receivers(dem, active2d, R, L)
    receiver = zeros(Int, length(dem))
    for c in R
        active2d[c] || continue
        z = dem[c]
        best_slope = 0.0
        best = 0
        for k in eachindex(_NB8)
            nc = c + _NB8[k]
            (nc in R && active2d[nc]) || continue
            drop = z - dem[nc]
            drop > 0 || continue
            slope = drop / _NB8_DIST[k]
            if slope > best_slope
                best_slope = slope
                best = L[nc]
            end
        end
        receiver[L[c]] = best
    end
    return receiver
end

# Priority-flood D8 (Geomorphometry): routes every pit up-and-over its spill
# point, so there are no terminal sinks. cellsize=(1,1) makes `_orient` the
# identity, so the LDD codes decode back to exact CartesianIndex receiver offsets.
function _priorityflood_receivers(dem, active2d, R, L)
    _, dirs = Geomorphometry.flowaccumulation(dem, .!active2d;
        method = D8(), cellsize = (1.0, 1.0))
    receiver = zeros(Int, length(dem))
    for c in R
        active2d[c] || continue
        off = CartesianIndex(dirs[c])
        off == CartesianIndex(0, 0) && continue
        rc = c + off
        (rc in R && active2d[rc]) && (receiver[L[c]] = L[rc])
    end
    return receiver
end

# Kahn topological order (sources first). Doubles as an acyclicity check: on a
# strictly-downhill forest every active cell is reachable, so the produced order
# must cover them all.
function _topological_order(receiver, indegree, active, ncells)
    order = Int[]
    sizehint!(order, count(active))
    remaining = copy(indegree)
    stack = Int[]
    for id in 1:ncells
        active[id] && indegree[id] == 0 && push!(stack, id)
    end
    while !isempty(stack)
        id = pop!(stack)
        push!(order, id)
        r = receiver[id]
        if r != 0
            remaining[r] -= 1
            remaining[r] == 0 && push!(stack, r)
        end
    end
    length(order) == count(active) || error(
        "build_flow_graph: flow graph is not acyclic (only strictly-downhill " *
        "routing is supported); $(count(active) - length(order)) cells are in cycles.")
    return order
end
