"""
    solve_max_dynamic_flow(net, s, t, T; method="FF", solver="ortools", solver_method="primal")

Solves the Maximum Dynamic Flow Problem (MDFP).

# Arguments
- `net`: The `NetworkInstance`.
- `s`, `t`: Source and sink nodes.
- `T`: The time horizon.
- `method`: `"FF"` (Ford-Fulkerson temporally repeated flow) or `"TE"` (Time-Expanded network).
- `solver`: `"HiGHS"` or `"Gurobi"`.
- `solver_method`: `"primal"`, `"dual"`, or `"network"`.
"""
function solve_max_dynamic_flow(net, s::Int, t::Int, T::Int; method::String="FF", solver::String="HiGHS", solver_method::String="primal")
    if method == "FF"
        session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)
        res = solve_mccp(session, T)
        
        if res.status == STATUS_OPTIMAL
            paths = decompose_flow_into_paths(net, s, t, res.flows)
            return MDFResult(STATUS_OPTIMAL, res.v, res.flows, paths)
        else
            return MDFResult(STATUS_INFEASIBLE, 0.0, Float64[], Tuple{Vector{Int},Float64}[])
        end
    elseif method == "TE"
        return solve_max_dynamic_flow_te_lp(net, s, t, T, solver, solver_method)
        # TODO: add the pure max flow implementation
    else
        error("Unknown method '$method'.")
    end
end
