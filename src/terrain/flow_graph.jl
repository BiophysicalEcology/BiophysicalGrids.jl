# flow_graph.jl
#
# A directed drainage graph built once from the run-grid DEM. Each active cell
# drains to a single downslope neighbour (D8), so the structure is a *forest*
# (out-degree ≤ 1) rooted at domain-edge / pit outlets — strictly acyclic for
# the downhill case we handle now.
#
# The graph is what lets a lateral process (overland flow, later cold-air
# drainage) stay compatible with the per-cell whole-time-series solve: cells are
# processed in dependency order (every cell after all its upslope contributors),
# so a downstream cell's inflow is the already-computed outflow time series of
# its upstream cells — read in as a forcing, no global time-stepping.
#
# Representation is three flat arrays (no per-node heap vectors, no external
# graph dependency): `receiver` (forward, single edge) and a CSR pair
# `upstream_ids` / `offsets` (inverse adjacency, contiguous per cell). This is
# the leanest form for the atomic dataflow scheduler in routing.jl, which
# indexes these arrays directly in its hot loop.

"""
    FlowGraph

Directed D8 drainage forest over the run grid. All ids are linear indices into
the **ascending (X→, Y→)** canonical matrix built from `elevation`; `dimindices`
maps each id back to the `(X(i), Y(j))` selector that indexes `terrain`,
`output`, and `build_inputs` (addressing is by dim name, so the storage order of
those arrays is irrelevant).

Fields:
- `receiver::Vector{Int}`     — downstream cell id per cell (`0` = outlet or inactive)
- `upstream_ids::Vector{Int}` — CSR store of contributor ids, grouped by receiver
- `offsets::Vector{Int}`      — CSR offsets; cell `id`'s contributors are
  `upstream_ids[offsets[id] : offsets[id+1]-1]` (length `ncells+1`)
- `indegree::Vector{Int}`     — number of upstream contributors (dependency count)
- `order::Vector{Int}`        — one valid topological order (sources first); used
  by tests/level analysis, not the dynamic scheduler
- `dimindices`                — cell id → `(X(i), Y(j))` selector
- `active::BitVector`         — active (unmasked, finite-elevation) cells
- `is_sink::BitVector`        — interior terminal cells (a local minimum fully
  surrounded by active cells; `receiver == 0`). Water routed here has nowhere
  to go and accumulates — an endorheic sink. Distinct from domain-edge outlets,
  which also have `receiver == 0` but let water leave the domain.
- `ncells::Int`
"""
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

# Contributor slice for cell `id` in the CSR store. Iterate as
# `for k in _upstream_range(g, id); up = g.upstream_ids[k]; ...; end`.
@inline _upstream_range(g::FlowGraph, id::Int) = g.offsets[id]:(g.offsets[id+1] - 1)

# Map from ascending-matrix position along one axis back to the raster's own
# index along that axis. Ascending axes map identically; descending axes flip.
@inline _asc_map(ascending::Bool, n::Int) = ascending ? (a -> a) : (a -> n - a + 1)

"""
    SinkHandling

How closed topographic depressions are treated when the flow graph is built —
the assumption about whether water reaching a local minimum retains or spills.
See [`TerminalSinks`](@ref) and [`SpillOver`](@ref).
"""
abstract type SinkHandling end

"""
    TerminalSinks(; pool = 1.0e12u"kg/m^2")

Closed depressions are **terminal endorheic sinks**: routed water accumulates in
the local minimum and can leave only *vertically* — through the column's own water
balance (evaporation, infiltration, deep drainage) — never laterally. Whether a
standing water body forms, and how deep, is an emergent result of that balance, not
imposed by the DEM. Flow directions are steepest-descent, so a cell with no
strictly-lower neighbour is terminal — which also makes undrained flats terminal
(a crude but conservative treatment). Default.

- `pool` — the `max_surface_pool` used for sink cells: large enough that routed
  inflow accumulates and is drawn down only by the column's evaporation and
  infiltration, instead of overflowing.
"""
@kwdef struct TerminalSinks <: SinkHandling
    pool::typeof(1.0u"kg/m^2") = 1.0e12u"kg/m^2"
end

"""
    SpillOver()

Closed depressions are assumed **already full and overflowing at their lowest
sill**, passing all inflow downslope with no net storage — a depressionless,
fully-integrated drainage network that reaches the domain edge everywhere. Assumes
an exorheic landscape (every basin drains out); nothing accumulates. Flow
directions come from a priority-flood that routes over every pit.
"""
struct SpillOver <: SinkHandling end

"""
    build_flow_graph(elevation::Raster, mask; method = D8(), sinks = TerminalSinks()) -> FlowGraph

Build the D8 drainage forest from the run-grid `elevation` raster (Unitful `m`)
and the active `mask` (a 2-D `Bool` raster on the same grid). Masked or
non-finite-elevation cells are excluded and never become a receiver.

`sinks::SinkHandling` controls what happens at closed depressions —
[`TerminalSinks`](@ref) (default, endorheic accumulation) or [`SpillOver`](@ref)
(priority-flood, drains through). Both are strictly downhill and acyclic, so the
scheduler and topological order are identical either way. Only single-receiver
`D8()` is supported.
"""
function build_flow_graph(elevation::Raster, mask; method = D8(), sinks::SinkHandling = TerminalSinks())
    method isa D8 || error(
        "build_flow_graph: only D8() routing is supported (single-receiver flow " *
        "with reconstructable weights); got $(method).")

    xlk = lookup(elevation, X)
    ylk = lookup(elevation, Y)
    nx = length(xlk)
    ny = length(ylk)
    xasc = nx < 2 || first(xlk) <= last(xlk)
    yasc = ny < 2 || first(ylk) <= last(ylk)
    ax = _asc_map(xasc, nx)
    ay = _asc_map(yasc, ny)

    # Ascending (X→, Y→) matrices: elevation (Float64, stripped) + active mask.
    dem = Matrix{Float64}(undef, nx, ny)
    active2d = falses(nx, ny)
    @inbounds for b in 1:ny, a in 1:nx
        e = elevation[X(ax(a)), Y(ay(b))]
        v = ustrip(u"m", e)
        on = isfinite(v) && Bool(mask[X(ax(a)), Y(ay(b))])
        active2d[a, b] = on
        dem[a, b] = on ? v : 0.0   # value under masked cells is never read by the flood
    end

    R = CartesianIndices(dem)
    L = LinearIndices(dem)
    ncells = length(dem)
    active = BitVector(vec(active2d))

    receiver = _receivers(sinks, dem, active2d, R, L)

    # Interior terminal cells (a local minimum fully surrounded by active,
    # in-grid cells) are endorheic sinks; terminal cells touching the domain
    # edge or a masked region are outlets (water leaves).
    is_sink = falses(ncells)
    @inbounds for c in R
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
    @inbounds for id in 1:ncells
        r = receiver[id]
        r == 0 && continue
        indegree[r] += 1
    end
    offsets = Vector{Int}(undef, ncells + 1)
    offsets[1] = 1
    @inbounds for id in 1:ncells
        offsets[id + 1] = offsets[id] + indegree[id]
    end
    upstream_ids = Vector{Int}(undef, offsets[end] - 1)
    cursor = copy(offsets)
    @inbounds for id in 1:ncells
        r = receiver[id]
        r == 0 && continue
        upstream_ids[cursor[r]] = id
        cursor[r] += 1
    end

    order = _topological_order(receiver, indegree, active, ncells)
    # Flat Vector indexed by linear cell id (column-major, matching `L`).
    dimindices = [(X(ax(c[1])), Y(ay(c[2]))) for c in vec(R)]

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
    @inbounds for c in R
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
    @inbounds for c in R
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
    @inbounds for id in 1:ncells
        active[id] && indegree[id] == 0 && push!(stack, id)
    end
    @inbounds while !isempty(stack)
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
