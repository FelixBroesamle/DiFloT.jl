const STATUS_OPTIMAL = 0
const STATUS_INFEASIBLE = 1

struct Edge
    tail::Int
    head::Int
    low::Int
    up::Int
    transit_time::Int
    parallel_index::Int
end

Edge(tail::Int, head::Int, low::Int, up::Int, transit_time::Int) =
    Edge(tail, head, low, up, transit_time, 0)

"""
    NetworkInstance

A static representation of a dynamic network instance.

# Fields
- `num_nodes`, `num_arcs`: sizes.
- `nodes`: node list.
- `edges`: arcs.
- `problem_type`: `:MIN`, `:MDFP`, or `:QFP`.
- `source`, `sink`: node ids (0 if unset).
- `target_value`: `T` for MDFP, `F` for QFP (0.0 if unset).
- `arc_indices`: quick lookup of arc positions.
"""
mutable struct NetworkInstance
    num_nodes::Int
    num_arcs::Int
    nodes::Vector{Int}
    edges::Vector{Edge}
    problem_type::Symbol
    source::Int
    sink::Int
    target_value::Float64
    arc_indices::Dict{Tuple{Int,Int},Vector{Int}}
end

"""
    NetworkInstance(
        start_nodes, end_nodes, capacities, transit_times;
        source=0, sink=1, target_value=1.0
    )

Array-based constructor.
"""
function NetworkInstance(
    start_nodes::Vector{Int}, 
    end_nodes::Vector{Int},
    capacities::Vector{Int}, 
    transit_times::Vector{Int};
    source::Int=0, 
    sink::Int=1, 
    target_value::Float64=1.0
)
    @assert length(start_nodes) == length(end_nodes) == length(capacities) == length(transit_times) "Array lengths must match."
    m = length(start_nodes)
    num_nodes = isempty(start_nodes) ? 0 : max(maximum(start_nodes), maximum(end_nodes))
    edges = Edge[]
    arc_indices = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in 1:m
        tail = start_nodes[i]
        head = end_nodes[i]
        k = (tail, head)
        slot = length(get(arc_indices, k, Int[]))
        push!(edges, Edge(tail, head, 0, capacities[i], transit_times[i], slot))
        push!(get!(arc_indices, k, Int[]), i)
    end
    return NetworkInstance(num_nodes, m, collect(1:num_nodes), edges,
                           :MIN, source, sink, target_value, arc_indices)
end

struct MCCPResult
    status::Int
    v::Float64
    g::Float64
    flows::Vector{Float64}
end

struct MDFResult
    status::Int
    value::Float64
    flows::Vector{Float64}
    paths::Vector{Tuple{Vector{Int},Float64}}
end

struct QFResult
    status::Int
    time_T::Float64
    static_flow::Vector{Float64}
end