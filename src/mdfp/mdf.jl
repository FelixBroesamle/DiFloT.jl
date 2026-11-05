"""
    solve_max_dynamic_flow(net, s, t, T; method="FF", solver="ortools", solver_method="primal")

Solves the Maximum Dynamic Flow Problem (MDFP).

# Arguments
- `net`: The `NetworkInstance`.
- `s`, `t`: Source and sink nodes.
- `T`: The time horizon.
- `method`: `"FF"` (Ford-Fulkerson temporally repeated flow) or `"TE"` (Time-Expanded network).
- `solver`: `"ortools"`, `"HiGHS"`, or `"Gurobi"`.
- `solver_method`: `"primal"`, `"dual"`, or `"network"`.
"""
function solve_max_dynamic_flow(net, s::Int, t::Int, T::Int; method::String="FF", solver::String="ortools", solver_method::String="primal")
    if method == "FF"
        if lowercase(solver) == "ortools"
            res = ORToolsBackend.solve_circulation(net, s, t, T)
        else
            session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)
            res = solve_mccp(session, T)
        end
        
        if res.status == STATUS_OPTIMAL
            paths = decompose_flow_into_paths(net, s, t, res.flows)
            return MDFResult(STATUS_OPTIMAL, res.v, res.flows, paths)
        else
            return MDFResult(STATUS_INFEASIBLE, 0.0, Float64[], Tuple{Vector{Int},Float64}[])
        end
    elseif method == "TE"
        if lowercase(solver) == "ortools"
            return solve_max_dynamic_flow_te(net, s, t, T)
        else
            return solve_max_dynamic_flow_te_lp(net, s, t, T, solver, solver_method)
        end
    else
        error("Unknown method '$method'.")
    end
end