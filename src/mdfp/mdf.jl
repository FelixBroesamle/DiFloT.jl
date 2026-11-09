"""
    solve_max_dynamic_flow(net, s, t, T; method="FF", solver="HiGHS", solver_method="primal")

Solves the Maximum Dynamic Flow Problem (MDFP).

# Arguments
- `net`: The `NetworkInstance`.
- `s`, `t`: Source and sink nodes.
- `T`: The time horizon.
- `method`: `"FF"` (Ford-Fulkerson) or `"TE"` (Time-Expanded network).
- `solver`: If `method="TE"`, passing `"HiGHS"` or `"Gurobi"` uses a Linear Program. Passing any other string (e.g., `"GraphsFlows"`) uses the highly optimized pure Julia Sparse Dinic algorithm.
- `solver_method`: `"primal"`, `"dual"`, or `"network"`.

# Returns
An `MDFResult` containing the `status`, `value` (maximum flow), and decomposed `paths`.
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
        #s_lower = lowercase(solver)
        #if s_lower == "highs" || s_lower == "gurobi"
        #    return solve_max_dynamic_flow_te_lp(net, s, t, T, solver, solver_method)
        #else
            # Defaults to GraphsFlows (Sparse Dinic)
        return solve_max_dynamic_flow_te(net, s, t, T)
        #end
    else
        error("Unknown method '$method'.")
    end
end