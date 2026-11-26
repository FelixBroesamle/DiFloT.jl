"""
    solve_qf_newton(net,s,t,F;solver="HiGHS",solver_method="primal", ...)

Newton-type search with the revised termination criterion |v(theta)-F|<1.
The real oracle uses the same model and cache as integer evaluations.
An invalid numerical step is reported rather than treated as convergence.
"""
function solve_qf_newton(net,s::Int,t::Int,F::Real;solver::String="HiGHS",
                         solver_method::String="primal",max_flow_solver::String="GraphsFlows",
                         max_iter::Int=1000,with_paths::Bool=false,verify::Bool=false)
    lowercase(solver) == "ortools" && throw(ArgumentError("Newton search requires a real-valued MCCP oracle."))
    F_int = _positive_demand(F)
    max_iter > 0 || throw(ArgumentError("newton_max_iter must be positive."))
    init = _initial_lower(net,s,t,F_int;solver=solver,solver_method=solver_method,max_flow_solver=max_flow_solver)
    isnothing(init) && return _infeasible_qf()
    (;session,T_l,r_l) = init
    _value(r_l) >= F_int && return _qf_result(session,T_l,r_l,F_int;with_paths=with_paths,verify=verify)
    theta = Float64(T_l)
    cur = r_l
    iter = 0
    while abs(cur.v-F_int) >= 1.0
        iter += 1
        iter <= max_iter || error("Newton search exceeded $max_iter iterations.")
        cur.g >= 1.0 || error("Invalid Newton subgradient $(cur.g).")
        next_theta = theta+(F_int-cur.v)/cur.g
        isfinite(next_theta) && next_theta >= T_l || error("Invalid Newton iterate $next_theta.")
        next_theta != theta || error("Newton search stagnated before meeting the residual criterion.")
        if iter > 1 && next_theta > theta
            error("Newton search lost monotonicity before convergence.")
        end
        theta = next_theta
        cur = solve_mccp(session,theta)
    end
    # Under the exact-oracle assumptions T* is floor(theta) or floor(theta)+1.
    T = floor(Int,theta)
    r = solve_mccp(session,T)
    predecessor_checked = _value(r) < F_int
    if predecessor_checked
        T += 1
        r = solve_mccp(session,T)
    end
    return _qf_result(session,T,r,F_int;with_paths=with_paths,verify=verify,predecessor_checked=predecessor_checked)
end
