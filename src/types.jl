"""Successful optimal solver status (LP status only for `LPCandidate`)."""
const STATUS_OPTIMAL = 0
"""Infeasible problem status."""
const STATUS_INFEASIBLE = 1

"""Directed arc with integer capacity, lower bound, transit time, and parallel-arc index."""
struct Arc
    tail::Int
    head::Int
    low::Int
    up::Int
    transit_time::Int
    parallel_index::Int
end

Arc(tail::Int, head::Int, low::Int, up::Int, transit_time::Int) = Arc(tail, head, low, up, transit_time, 0)

"""
    NetworkInstance

A static representation of a dynamic network instance.

# Fields
- `num_nodes`, `num_arcs`: Network size.
- `nodes`: Node IDs in `1:num_nodes`.
- `arcs`: Directed arcs in their original order.
- `problem_type`: `:MIN`, `:MDFP`, or `:QFP`.
- `source`, `sink`: Terminal IDs (0 if unset).
- `target_value`: Time-horizon T or demand F.
- `arc_indices`: Global arc positions for each directed pair `(tail, head)`.

Parallel arcs are distinguished by a zero-based `parallel_index`.
Global arc positions are one-based.
"""
mutable struct NetworkInstance
    num_nodes::Int
    num_arcs::Int
    nodes::Vector{Int}
    arcs::Vector{Arc}
    problem_type::Symbol
    source::Int
    sink::Int
    target_value::Float64
    arc_indices::Dict{Tuple{Int,Int},Vector{Int}}
end

"""
    NetworkInstance(
        start_nodes, end_nodes, capacities, transit_times;
        source=1, sink=2, target_value=1.0
    )

Construct a network from arrays. Lower bounds are zero.
"""
function NetworkInstance(
    start_nodes::Vector{Int}, end_nodes::Vector{Int}, capacities::Vector{Int}, transit_times::Vector{Int};
    source::Int=1, sink::Int=2, target_value::Real=1.0
)
    length(start_nodes) == length(end_nodes) == length(capacities) == length(transit_times) ||
        throw(ArgumentError("Array lengths must match."))
    
    all(>(0), start_nodes) && all(>(0), end_nodes) ||
        throw(ArgumentError("Node IDs must be positive."))
    
    source >= 0 && sink >= 0 || throw(ArgumentError("Invalid terminal ID."))
    all(>=(0), capacities) || throw(ArgumentError("Negative capacity."))
    all(>=(0), transit_times) || throw(ArgumentError("Negative transit time."))
    
    m = length(start_nodes)
    n = max(source, sink, maximum(start_nodes; init=0), maximum(end_nodes; init=0))
    
    arcs = Arc[]
    arc_indices = Dict{Tuple{Int,Int},Vector{Int}}()
    
    for i in 1:m
        tail, head = start_nodes[i], end_nodes[i]
        indices = get!(arc_indices, (tail, head), Int[])
    
        push!(arcs, Arc(tail, head, 0, capacities[i], transit_times[i], length(indices)))
        push!(indices, i)
    end
    
    return NetworkInstance(
        n, m, collect(1:n), arcs, :MIN, source, sink, _target_float(target_value), arc_indices
    )
end

const PathFlow = Tuple{Vector{Int},Float64}

struct MCCPResult
    status::Int
    v::Float64
    g::Float64
    flows::Vector{Float64}
end

"""
    MDFResult

Result of the Maximum Dynamic Flow Problem.
- `status`: Solution status.
- `value`: Maximum dynamic flow value.
- `flows`: Static arc flows (FF only).
- `paths`: Optional path decomposition, using global arc indices.
"""
struct MDFResult
    status::Int
    value::Float64
    flows::Vector{Float64}
    paths::Vector{PathFlow}
end

"""
    QFResult

Result of the Quickest Flow Problem.
- `status`: Solution status.
- `time_T`: Minimum time-horizon.
- `static_flow`: Corresponding optimal static flow.
- `paths`: Optional path decomposition.
- `mccp_calls`: Number of actual MCCP oracle calls, including initialization.
- `lp_calls`: Number of LP oracle calls.
- `safeguard_activations`: Number of midpoint fallbacks (INT-SG).
"""
struct QFResult
    status::Int
    time_T::Float64
    static_flow::Vector{Float64}
    paths::Vector{PathFlow}
    mccp_calls::Int
    lp_calls::Int
    safeguard_activations::Int
end

QFResult(status, time_T, static_flow) = QFResult(status, time_T, static_flow, PathFlow[], 0, 0, 0)

"""
    LPCandidate

Numerical Lin--Jaillet LP solution and rounded discrete time-horizon.
`status` refers to the LP, not the optimality of the discrete candidate.
"""
struct LPCandidate
    status::Int
    tau_tilde::Float64
    time_T::Int
    lp_calls::Int
end
