"""Quickest Flow via improved binary search (secant and supporting lines)."""
function solve_qf_ibs(
    net, s::Int, t::Int, F::Real; 
    solver::String="HiGHS", solver_method::String="primal", max_flow_solver::String="GraphsFlows",
    with_paths::Bool=false, verify::Bool=false
)
    F_int = _positive_demand(F)

    # Initialize the lower bound and its MCCP solution.
    init = _initial_lower(
        net, s, t, F_int;
        solver=solver, solver_method=solver_method, max_flow_solver=max_flow_solver
    )
    
    isnothing(init) && return _infeasible_qf()
    (; session, T_l, r_l) = init

    # The lower bound is optimal if it already satisfies the demand.
    if _value(r_l) >= F_int
        return _qf_result(
            session, T_l, r_l, F_int;
            with_paths=with_paths, verify=verify
        )
    end

    # Establish the bracket v(T_l) < F <= v(T_u).
    T_u, r_u = _initial_upper(session, T_l, r_l, F_int)

    while T_l + 1 < T_u
        v_l, g_l = _value(r_l), _slope(r_l)
        v_u, g_u = _value(r_u), _slope(r_u)
        g_l >= 1 && g_u >= g_l && v_l < F_int <= v_u || error("Invalid IBS bracket or slopes.")

        # Lower bound: ceiling of the secant root (exact arithmetic).
        alpha = T_l + Int(cld(Int128(T_u - T_l) * (F_int - v_l), Int128(v_u - v_l)))

        # Upper bound: the smaller ceiling of the supporting-line roots.
        sigma_l = T_l + cld(F_int - v_l, g_l)
        sigma_u = T_u + cld(F_int - v_u, g_u)
        beta = min(sigma_l, sigma_u)
        T_l + 1 <= alpha <= beta <= T_u || error("Inconsistent IBS bounds ($alpha,$beta).")

        # The tightened bounds may already identify the optimum.
        if alpha == beta
            r = alpha == T_u ? r_u : solve_mccp(session, alpha)
            return _qf_result(
                session, alpha, r, F_int;
                with_paths=with_paths, verify=verify
            )
        end

        # Bisect the tightened interval.
        T_m = _midpoint(alpha, beta)
        r_m = solve_mccp(session, T_m)
        
        if _value(r_m) >= F_int
            T_u, r_u = T_m, r_m
        else
            T_l, r_l = T_m, r_m
        end
    end
    
    return _qf_result(
        session, T_u, r_u, F_int;              
        with_paths=with_paths, verify=verify, predecessor_checked=true
    )
end
