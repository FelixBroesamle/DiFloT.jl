"""Return the unique arc key (tail, head, parallel index)."""
arc_key(e::Arc) = (e.tail, e.head, e.parallel_index)

"""Return incoming and outgoing arc keys for each node."""
function incidence_mapping(net::NetworkInstance)
    incoming = [Tuple{Int,Int,Int}[] for _ in 1:net.num_nodes]
    outgoing = [Tuple{Int,Int,Int}[] for _ in 1:net.num_nodes]
    for e in net.arcs
        push!(outgoing[e.tail], arc_key(e))
        push!(incoming[e.head], arc_key(e))
    end
    return incoming, outgoing
end

"""Return incoming and outgoing neighboring nodes for each node."""
function adjacency_mapping(net::NetworkInstance)
    incoming = [Int[] for _ in 1:net.num_nodes]
    outgoing = [Int[] for _ in 1:net.num_nodes]
    for e in net.arcs
        push!(outgoing[e.tail], e.head)
        push!(incoming[e.head], e.tail)
    end
    return incoming, outgoing
end

"""Return incoming and outgoing arc positions for each node."""
function positional_incidence(net::NetworkInstance)
    incoming = [Int[] for _ in 1:net.num_nodes]
    outgoing = [Int[] for _ in 1:net.num_nodes]
    for (i,e) in enumerate(net.arcs)
        push!(outgoing[e.tail], i)
        push!(incoming[e.head], i)
    end
    return incoming, outgoing
end
