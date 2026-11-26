# The public result types retain Float64 fields. Restrict integer quantities
# to a range in which all integer values can be represented exactly.
const MAX_EXACT_INT = Int64(2)^53 - 1
const FLOW_TOL = 1e-5

function _exact_int(x::Real, name::AbstractString)
    isfinite(x) && isinteger(x) && abs(x) <= MAX_EXACT_INT ||
        throw(ArgumentError("$name must be an integer in the supported exact Float64 range; got $x."))
    return Int(x)
end

function _target_float(x::Real)
    isfinite(x) && abs(x) <= MAX_EXACT_INT || throw(ArgumentError("Invalid target value."))
    return Float64(x)
end

function _positive_demand(F::Real)
    f = _exact_int(F, "F")
    f > 0 || throw(ArgumentError("F must be positive."))
    return f
end

function _validate_terminals(net::NetworkInstance, s::Int, t::Int)
    1 <= s <= net.num_nodes && 1 <= t <= net.num_nodes && s != t ||
        throw(ArgumentError("Source and sink must be distinct valid node IDs."))
    return nothing
end

"""Validate the network assumptions used in the paper."""
function _validate_network(net::NetworkInstance, s::Int, t::Int)
    _validate_terminals(net, s, t)
    net.num_arcs == length(net.arcs) || throw(ArgumentError("Inconsistent arc count."))
    total_capacity = Int128(0)
    total_cost = Int128(0)
    for e in net.arcs
        1 <= e.tail <= net.num_nodes && 1 <= e.head <= net.num_nodes ||
            throw(ArgumentError("Invalid arc endpoint."))
        e.low == 0 || throw(ArgumentError("Dynamic flow methods require zero lower bounds."))
        0 <= e.up <= MAX_EXACT_INT && 0 <= e.transit_time <= MAX_EXACT_INT ||
            throw(ArgumentError("Invalid capacity or transit time."))
        total_capacity += e.up
        total_cost += Int128(e.up) * e.transit_time
        total_capacity <= MAX_EXACT_INT && total_cost <= MAX_EXACT_INT ||
            throw(ArgumentError("Network exceeds the supported exact-integer range."))
    end
    return nothing
end

function _check_horizon(net::NetworkInstance, s::Int, t::Int, T::Int)
    0 <= T < MAX_EXACT_INT || throw(ArgumentError("Invalid time-horizon $T."))
    Int128(T + 1) * c_rate_bound(net, s, t) <= MAX_EXACT_INT ||
        throw(ArgumentError("Horizon/capacity product exceeds the supported exact-integer range."))
    return nothing
end

function _rounded_flow(x::Real, label::AbstractString; tol::Real=FLOW_TOL)
    isfinite(x) && abs(x) <= MAX_EXACT_INT || error("Invalid $label: $x.")
    y = round(Int, x)
    abs(x-y) <= tol || error("Nonintegral $label: $x. Check solver accuracy and algorithm.")
    return y
end

function _require_optimal(model::JuMP.Model, name::AbstractString)
    st = termination_status(model)
    st == MOI.OPTIMAL && JuMP.primal_status(model) == MOI.FEASIBLE_POINT ||
        error("$name failed: termination status $st, primal status $(JuMP.primal_status(model)).")
    return nothing
end
