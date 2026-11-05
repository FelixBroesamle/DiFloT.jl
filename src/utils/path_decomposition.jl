"""
    decompose_flow_into_paths(net, s, t, flow) -> Vector{Tuple{Vector{Int}, Float64}}

Decomposes a static flow vector into a set of s-t paths and their flow values.
"""
function decompose_flow_into_paths(net::NetworkInstance, s::Int, t::Int, flow::Vector{Float64})
    f = copy(flow)
    paths = Tuple{Vector{Int},Float64}[]
    _, out_edges = positional_incidence(net)

    while true
        visited     = fill(false, net.num_nodes)
        parent_edge = fill(-1, net.num_nodes)
        stack       = [s]
        visited[s]  = true
        reached_t   = false

        while !isempty(stack) && !reached_t
            u = pop!(stack)
            for ei in out_edges[u]
                if f[ei] > 1e-6
                    v = net.edges[ei].head
                    if !visited[v]
                        visited[v]     = true
                        parent_edge[v] = ei
                        push!(stack, v)
                        if v == t
                            reached_t = true
                            break
                        end
                    end
                end
            end
        end
        reached_t || break

        curr, bottleneck, path_edges = t, Inf, Int[]
        while curr != s
            ei = parent_edge[curr]
            pushfirst!(path_edges, ei)
            bottleneck = min(bottleneck, f[ei])
            curr       = net.edges[ei].tail
        end
        push!(paths, (path_edges, bottleneck))
        for ei in path_edges
            f[ei] -= bottleneck
        end
    end
    return paths
end