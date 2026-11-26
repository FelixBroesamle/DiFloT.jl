_infeasible_qf() = QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
_value(r::MCCPResult) = Int(r.v)  # only for integer-horizon evaluations
_slope(r::MCCPResult) = Int(r.g)
_midpoint(T_l::Int, T_u::Int) = T_l + fld(T_u - T_l, 2)

"""Initial lower bound, with MCCP(0) evaluated only if T_0=0."""
function _initial_lower(net, s, t, F::Int; solver, solver_method, max_flow_solver)
    # Shortest-path length and maximum static flow rate.
    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return nothing
    T0 = Int(T0f)
    
    gmax = compute_max_static_flow(net, s, t; solver=max_flow_solver).value
    gmax == 0 && return nothing

    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)

    # v(0) is zero unless a zero-transit-time path exists.
    v0 = T0 == 0 ? _value(solve_mccp(session, 0)) : 0

    # Combine the shortest-path and flow-rate lower bounds.
    T_l = max(T0, cld(F - v0, gmax))
    r_l = solve_mccp(session, T_l)
    
    return (session=session, T_l=T_l, r_l=r_l, gmax=gmax)
end

"""Initial feasible upper bound via a supporting line."""
function _initial_upper(session, T_l::Int, r_l::MCCPResult, F::Int)
    g_l = _slope(r_l)
    g_l >= 1 || error("Nonpositive slope at initial lower bound.")

    # Round up the intersection of the supporting line with F.
    T_u = T_l + cld(F - _value(r_l), g_l)
    r_u = solve_mccp(session, T_u)
    _value(r_l) < F <= _value(r_u) || error("Invalid initial QFP bracket.")
    
    return T_u, r_u
end

"""Build the result; optional additional predecessor verification is counted."""
function _qf_result(
    session, T::Int, r::MCCPResult, F::Int;
    with_paths::Bool=false, verify::Bool=false, lp_calls::Int=0, safeguard_activations::Int=0, predecessor_checked::Bool=false
)
    _value(r) >= F || error("Returned horizon $T does not satisfy demand $F.")

    # Optional predecessor check; count any additional oracle solve.
    if verify && T > 0 && !predecessor_checked
        _value(solve_mccp(session, T - 1)) < F ||
            error("Returned horizon $T is not minimal according to the MCCP oracle.")
    end
    
    paths = with_paths ? decompose_flow_into_paths(session.net, session.s, session.t, r.flows; T=T) : PathFlow[]
    
    return QFResult(
        STATUS_OPTIMAL, Float64(T), r.flows, paths, session.calls, lp_calls, safeguard_activations
    )
end
