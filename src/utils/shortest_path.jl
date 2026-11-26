"""
    compute_shortest_path(net, s, t) -> Float64

Shortest transit time over positive-capacity arcs; `Inf` if unreachable.
The graph topology is built independently of the weight matrix to preserve
zero-transit-time arcs.
"""
function compute_shortest_path(net::NetworkInstance, s::Int, t::Int)
    if s == t
        1 <= s <= net.num_nodes || throw(ArgumentError("Invalid node ID."))
        return 0.0
    end
    _validate_network(net,s,t)
    n = net.num_nodes
    g = Graphs.SimpleDiGraph(n)
    weights = Dict{Tuple{Int,Int},Int}()
    for e in net.arcs
        e.up > 0 && e.tail != e.head || continue
        k = (e.tail,e.head)
        weights[k] = min(get(weights,k,typemax(Int)),e.transit_time)
        add_edge!(g,e.tail,e.head)
    end
    I, J, V = Int[], Int[], Float64[]
    for ((u,v),w) in weights
        push!(I,u); push!(J,v); push!(V,Float64(w))
    end
    dist = sparse(I,J,V,n,n)
    d = Graphs.dijkstra_shortest_paths(g,s,dist).dists[t]
    return isfinite(d) && d != typemax(Float64) ? Float64(d) : Inf
end
