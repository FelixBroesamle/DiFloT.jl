using SimpleWeightedGraphs

"""
    compute_shortest_path(net::NetworkInstance, s::Int, t::Int) -> Float64

Computes the shortest path from `s` to `t` with respect to `transit_time`.
Returns `Inf` if no path exists.
"""
function compute_shortest_path(net::NetworkInstance, s::Int, t::Int)
    if s == t; return 0.0; end
    
    n = net.num_nodes
    
    # 1. Find the minimum transit time for parallel arcs
    min_weights = Dict{Tuple{Int,Int}, Float64}()
    for e in net.arcs
        k = (e.tail, e.head)
        w = Float64(e.transit_time)
        if !haskey(min_weights, k) || w < min_weights[k]
            min_weights[k] = w
        end
    end

    # 2. Extract arrays and build the sparse matrix directly
    sources = [k[1] for k in keys(min_weights)]
    destinations = [k[2] for k in keys(min_weights)]
    weights = [min_weights[k] for k in keys(min_weights)]

    # Build sparse matrix explicitly with size n x n to handle isolated nodes
    mat = sparse(sources, destinations, weights, n, n)

    # 3. Construct the graph in one shot from the sparse matrix
    g = SimpleWeightedDiGraph(mat)

    # 4. Run Dijkstra
    state = Graphs.dijkstra_shortest_paths(g, s)
    d = state.dists[t]
    
    return isfinite(d) ? Float64(d) : Inf
end