using GraphsFlows

"""
    compute_max_static_flow(net, s, t)

Computes the maximum static flow value `g_max` using the pure Julia `GraphsFlows.jl` backend
with Sparse Matrices and Dinic's Algorithm for optimal performance.

# Arguments
- `net`: The `NetworkInstance` containing the network topology.
- `s`: The source node ID.
- `t`: The sink node ID.

# Returns
A `NamedTuple` with fields `status` (`STATUS_OPTIMAL` or `STATUS_INFEASIBLE`) and `value` (the maximum flow).
"""
function compute_max_static_flow(net::NetworkInstance, s::Int, t::Int)
    n = net.num_nodes
    
    # Build the sparse capacity matrix directly from the arcs
    I = Int[]
    J = Int[]
    V = Float64[]
    
    for e in net.arcs
        if e.tail != e.head && e.up > 0
            push!(I, e.tail)
            push!(J, e.head)
            push!(V, Float64(e.up))
        end
    end
    
    # Sparse automatically sums capacities for parallel arcs!
    capacity_matrix = sparse(I, J, V, n, n)
    flow_graph = Graphs.SimpleDiGraph(capacity_matrix)
    
    # Run Dinic's Algorithm on the sparse matrix
    f_val, _ = GraphsFlows.maximum_flow(flow_graph, s, t, capacity_matrix, algorithm=GraphsFlows.PushRelabelAlgorithm())
    return (status = STATUS_OPTIMAL, value = f_val)
end