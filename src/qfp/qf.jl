"""
    solve_quickest_flow(
        net, s, t, F; 
        method="INT-SG", solver="HiGHS", solver_method="primal", max_flow_solver="GraphsFlows",
        max_iter=100, tol=1e-7, with_paths=false, verify=false, newton_max_iter=1000, mccp_solver=solver,
        mccp_solver_method="primal"
    )

Solves the Quickest Flow Problem (QFP).

# Arguments
- `net`: The `NetworkInstance`.
- `s`, `t`: Source and sink nodes.
- `F`: Positive integer demand.
- `method`: `"BS"`, `"IBS"`, `"INT"`, `"INT-SG"`, `"NEWT"`, or `"LP"` (also `"LP-check"`).
- `solver`: `"HiGHS"`, `"Gurobi"`, or `"ORTools"` (BS, IBS, INT, INT-SG only).
- `solver_method`: `"primal"` or `"dual"`.
- `max_flow_solver`: `"GraphsFlows"` or `"ORTools"` for initialization.
- `max_iter`, `tol`: Parameters for interpolation root-finding.
- `with_paths`: Include the path decomposition.
- `verify`: Additionally check the predecessor horizon, if not already checked.
- `newton_max_iter`: Safety cap for Newton search.
- `mccp_solver`, `mccp_solver_method`: MCCP backend for LP recovery.

# Returns
A `QFResult` with the solution status, minimum time-horizon, corresponding
static flow, and number of actual MCCP oracle calls. Additional verification
calls are counted. The LP approach assumes objective error < 1, as in the
paper. Solver failures raise errors; infeasible networks return `STATUS_INFEASIBLE`.
"""
function solve_quickest_flow(
    net, s::Int, t::Int, F::Real;
    method::String="INT-SG", solver::String="HiGHS", solver_method::String="primal", max_flow_solver::String="GraphsFlows",
    max_iter::Int=100, tol::Float64=1e-7, with_paths::Bool=false, verify::Bool=false,
    newton_max_iter::Int=1000, mccp_solver::String=solver, mccp_solver_method::String="primal"
)
    F_int = _positive_demand(F)
    _validate_network(net, s, t)
    method in ("BS", "IBS", "INT", "INT-SG", "NEWT", "LP", "LP-check") ||
        throw(ArgumentError("Unknown method '$method'."))

    # OR-Tools MCCP costs must be integral; the LP requires an LP solver.
    if lowercase(solver) == "ortools" && method in ("NEWT", "LP", "LP-check")
        throw(ArgumentError("Method '$method' requires HiGHS or Gurobi."))
    end

    # Shared options for the selected solution method.
    common = (; solver, solver_method, max_flow_solver, with_paths, verify)
    
    if method == "BS"
        return solve_qf_bs(net, s, t, F_int; common...)
    elseif method == "IBS"
        return solve_qf_ibs(net, s, t, F_int; common...)
    elseif method == "INT"
        return solve_qf_intp(net, s, t, F_int; common..., max_iter=max_iter, tol=tol)
    elseif method == "INT-SG"
        return solve_qf_intsg(net, s, t, F_int; common..., max_iter=max_iter, tol=tol)
    elseif method == "NEWT"
        return solve_qf_newton(net, s, t, F_int; common..., max_iter=newton_max_iter)
    end

    return solve_qf_lp(
        net, s, t, F_int; 
        common...,
        mccp_solver=mccp_solver, mccp_solver_method=mccp_solver_method
    )
end
