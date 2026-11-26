function _max_flow_value(I::Vector{Int}, J::Vector{Int}, V::Vector{Int},
                         n::Int, s::Int, t::Int, solver::String)
    if lowercase(solver) == "ortools"
        return _ortools_request(1,I,J,V,n,s,t).value
    elseif lowercase(solver) == "graphsflows"
        # Sparse construction sums capacities of parallel arcs.
        capacity_matrix = sparse(I,J,V,n,n)
        flow_graph = Graphs.SimpleDiGraph(capacity_matrix)
        f_val, _ = GraphsFlows.maximum_flow(flow_graph,s,t,capacity_matrix;
                                          algorithm=GraphsFlows.PushRelabelAlgorithm())
        return _rounded_flow(f_val, "maximum-flow value")
    end
    throw(ArgumentError("Unknown max-flow solver '$solver'."))
end

"""
    compute_max_static_flow(net, s, t; solver="GraphsFlows")

Computes the integral maximum static flow value using push-relabel.

# Arguments
- `net`: The `NetworkInstance`.
- `s`, `t`: Source and sink nodes.
- `solver`: `"GraphsFlows"` or `"ORTools"`.

# Returns
A named tuple with fields `status` and `value`.
A value of zero is a valid optimum, not a solver failure.
"""
function compute_max_static_flow(net::NetworkInstance, s::Int, t::Int; solver::String="GraphsFlows")
    _validate_network(net,s,t)
    I, J, V = Int[], Int[], Int[]
    for e in net.arcs
        if e.tail != e.head && e.up > 0
            push!(I,e.tail); push!(J,e.head); push!(V,e.up)
        end
    end
    (isempty(V) || c_rate_bound(net,s,t) == 0) && return (status=STATUS_OPTIMAL,value=0)
    f = _max_flow_value(I,J,V,net.num_nodes,s,t,solver)
    0 <= f <= c_rate_bound(net,s,t) || error("Invalid maximum-flow value $f.")
    return (status=STATUS_OPTIMAL,value=f)
end
