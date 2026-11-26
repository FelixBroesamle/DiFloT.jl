"""
    decompose_flow_into_paths(net, s, t, flow; T=nothing)

Decompose the positive net flow into simple s-t paths, represented by arc
positions. Discard residual circulations and, if T is given, paths whose
transit time exceeds T. The algorithm extracts at most m paths.
"""
function decompose_flow_into_paths(net::NetworkInstance,s::Int,t::Int,flow::Vector{Float64};
                                   T::Union{Nothing,Int}=nothing, tol::Float64=1e-5)
    _validate_terminals(net,s,t)
    length(flow) == net.num_arcs || throw(ArgumentError("Invalid flow length."))
    
    # The MCCP backend supplies integral flows; also accept fractional flows.
    isfinite(tol) && tol >= 0 || throw(ArgumentError("Invalid flow tolerance."))
    all(isfinite,flow) || error("Nonfinite static flow.")
    
    f = copy(flow)
    
    # Compute the net static flow value at the source.
    remaining = 0.0
    for (i, e) in enumerate(net.arcs)
        if e.tail == s
            remaining += f[i]
        end
        if e.head == s
            remaining -= f[i]
        end
    end

    paths = PathFlow[]
    _, out_edges = positional_incidence(net)
    visited = falses(net.num_nodes)
    parent_edge = zeros(Int,net.num_nodes)
    stack = Int[]
    
    while remaining > tol
        fill!(visited,false); empty!(stack)
        push!(stack,s); visited[s] = true
    
        while !isempty(stack) && !visited[t]
            u = pop!(stack)
    
            for i in out_edges[u]
                v = net.arcs[i].head
                if f[i] > 0.0 && !visited[v]
                    visited[v] = true; parent_edge[v] = i
                    push!(stack,v)
                    v == t && break
                end
            end
        end

        visited[t] || error("Positive residual value but no s-t path.")
        curr, delta, path = t, remaining, Int[]
        
        transit = 0
        while curr != s
            i = parent_edge[curr]
            push!(path,i)
            delta = min(delta,f[i])
            transit += net.arcs[i].transit_time
            curr = net.arcs[i].tail
        end
        
        reverse!(path)
        if isnothing(T) || transit <= T
            push!(paths,(path,Float64(delta)))
        end
        
        for i in path
            f[i] -= delta
        end
        
        remaining -= delta
    end
    return paths
end
