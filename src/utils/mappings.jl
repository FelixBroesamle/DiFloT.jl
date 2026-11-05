"""
    arc_key(e::Edge) -> Tuple{Int,Int,Int}

Returns the unique `(tail, head, parallel_index)` tuple for an edge.
"""
arc_key(e::Edge) = (e.tail, e.head, e.parallel_index)

"""
    incidence_mapping(net::NetworkInstance) -> (incoming, outgoing)

Returns dictionaries mapping each node to its incoming and outgoing arc keys.
"""
function incidence_mapping(net::NetworkInstance)
    incoming = [Tuple{Int,Int,Int}[] for _ in 1:net.num_nodes]
    outgoing = [Tuple{Int,Int,Int}[] for _ in 1:net.num_nodes]
    for e in net.edges
        k = arc_key(e)
        push!(outgoing[e.tail], k)
        push!(incoming[e.head], k)
    end
    return incoming, outgoing
end

"""
    adjacency_mapping(net::NetworkInstance) -> (incoming, outgoing)

Returns arrays mapping each node to its neighboring nodes (tails and heads).
"""
function adjacency_mapping(net::NetworkInstance)
    incoming = [Int[] for _ in 1:net.num_nodes]
    outgoing = [Int[] for _ in 1:net.num_nodes]
    for e in net.edges
        push!(outgoing[e.tail], e.head)
        push!(incoming[e.head], e.tail)
    end
    return incoming, outgoing
end

"""
    positional_incidence(net::NetworkInstance) -> (incoming, outgoing)

Returns arrays mapping each node to the positions (indices) of its incident arcs in `net.edges`.
"""
function positional_incidence(net::NetworkInstance)
    incoming = [Int[] for _ in 1:net.num_nodes]
    outgoing = [Int[] for _ in 1:net.num_nodes]
    for (i, e) in enumerate(net.edges)
        push!(outgoing[e.tail], i)
        push!(incoming[e.head], i)
    end
    return incoming, outgoing
end