"""Upper bound on the static flow rate (the smaller terminal cut)."""
function c_rate_bound(net::NetworkInstance, s::Int, t::Int)
    out_s = sum(Int128(e.up) for e in net.arcs if e.tail == s && e.head != s; init=Int128(0))
    in_t = sum(Int128(e.up) for e in net.arcs if e.head == t && e.tail != t; init=Int128(0))
    
    return Int(min(out_s, in_t))
end
