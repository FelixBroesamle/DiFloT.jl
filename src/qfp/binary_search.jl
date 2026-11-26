"""Quickest Flow via standard binary search."""
function solve_qf_bs(
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

    # Bisect until the infeasible and feasible horizons are consecutive.
    while T_l + 1 < T_u
        T_m = _midpoint(T_l, T_u)
        r_m = solve_mccp(session, T_m)

        if _value(r_m) >= F_int
            T_u, r_u = T_m, r_m
        else
            T_l, r_l = T_m, r_m
        end
    end

    # The feasible endpoint is the minimum time-horizon.
    return _qf_result(
        session, T_u, r_u, F_int;
        with_paths=with_paths, verify=verify, predecessor_checked=true
    )
end