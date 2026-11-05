"""
    c_rate_bound(net::NetworkInstance, s::Int, t::Int) -> Int

Computes the upper bound on the static max flow `g_max`, which is the minimum of the total capacity
leaving the source and the total capacity entering the sink.
"""
function c_rate_bound(net::NetworkInstance, s::Int, t::Int)
    out_s = sum(e.up for e in net.edges if e.tail == s; init = 0)
    in_t  = sum(e.up for e in net.edges if e.head == t; init = 0)
    return min(out_s, in_t)
end