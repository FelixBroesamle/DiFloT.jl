"""
    compute_shortest_path(net::NetworkInstance, s::Int, t::Int) -> Float64

Computes the shortest path from `s` to `t` with respect to `transit_time`.
Returns `Inf` if no path exists.
"""
function compute_shortest_path(net::NetworkInstance, s::Int, t::Int)
    if s == t; return 0.0; end
    g = Graphs.SimpleDiGraph(net.num_nodes)
    distmx = fill(Inf, net.num_nodes, net.num_nodes)
    for e in net.edges
        Graphs.add_edge!(g, e.tail, e.head)
        c = Float64(e.transit_time)
        if c < distmx[e.tail, e.head]
            distmx[e.tail, e.head] = c
        end
    end
    state = Graphs.dijkstra_shortest_paths(g, s, distmx)
    d = state.dists[t]
    return isfinite(d) ? Float64(d) : Inf
end