"""Quickest Flow via standard binary search."""
function solve_qf_bs(net,s::Int,t::Int,F::Real;solver::String="HiGHS",
                     solver_method::String="primal",max_flow_solver::String="GraphsFlows",
                     with_paths::Bool=false,verify::Bool=false)
    F_int = _positive_demand(F)
    init = _initial_lower(net,s,t,F_int;solver=solver,solver_method=solver_method,max_flow_solver=max_flow_solver)
    isnothing(init) && return _infeasible_qf()
    (;session,T_l,r_l) = init
    _value(r_l) >= F_int && return _qf_result(session,T_l,r_l,F_int;with_paths=with_paths,verify=verify)
    T_u,r_u = _initial_upper(session,T_l,r_l,F_int)
    T_u,r_u = _binary_search(session,T_l,r_l,T_u,r_u,F_int)
    return _qf_result(session,T_u,r_u,F_int;with_paths=with_paths,verify=verify,predecessor_checked=true)
end
