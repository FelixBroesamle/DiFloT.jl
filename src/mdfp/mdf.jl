"""
    solve_max_dynamic_flow(net, s, t, T; method="FF", solver="HiGHS",
                           solver_method="primal", max_flow_solver="GraphsFlows",
                           with_paths=true)

Solves the Maximum Dynamic Flow Problem (MDFP).

# Arguments
- `net`: The `NetworkInstance`.
- `s`, `t`: Source and sink nodes.
- `T`: Nonnegative integer time-horizon.
- `method`: `"FF"` (Ford--Fulkerson) or `"TE"` (time-expanded network).
- `solver`: MCCP backend (`"HiGHS"`, `"ORTools"`, or `"Gurobi"`).
- `solver_method`: `"primal"` or `"dual"` for LP-based backends.
- `max_flow_solver`: Maximum-flow backend (`"GraphsFlows"` or `"ORTools"`).
- `with_paths`: Include the path decomposition for FF.

# Returns
An `MDFResult` containing the maximum flow value and, for FF, static arc
flows and optional paths. TE returns the value only.
"""
function solve_max_dynamic_flow(net,s::Int,t::Int,T::Int;
                                method::String="FF",solver::String="HiGHS",
                                solver_method::String="primal",
                                max_flow_solver::String="GraphsFlows",
                                with_paths::Bool=true)
    _validate_network(net,s,t)
    _check_horizon(net,s,t,T)
    method in ("FF","TE") || throw(ArgumentError("Unknown method '$method'."))
    T0 = compute_shortest_path(net,s,t)
    if isinf(T0) || T0 > T
        flows = method == "FF" ? zeros(Float64,net.num_arcs) : Float64[]
        return MDFResult(STATUS_OPTIMAL,0.0,flows,PathFlow[])
    end
    if method == "FF"
        if lowercase(solver) == "ortools"
            res = _ortools_mccp(net,s,t,T)
        else
            session = MCCPSession(net,s,t;solver=solver,solver_method=solver_method)
            res = solve_mccp(session,T)
        end
        paths = with_paths ? decompose_flow_into_paths(net,s,t,res.flows;T=T) : PathFlow[]
        return MDFResult(STATUS_OPTIMAL,res.v,res.flows,paths)
    end
    # Preserve TE's original combinatorial default; solver can also select
    # a max-flow backend explicitly for convenience.
    backend = lowercase(solver) in ("graphsflows","ortools") ? solver : max_flow_solver
    return solve_max_dynamic_flow_te(net,s,t,T;solver=backend)
end

"""Alias for [`solve_max_dynamic_flow`](@ref)."""
solve_maximum_dynamic_flow(args...;kwargs...) = solve_max_dynamic_flow(args...;kwargs...)
