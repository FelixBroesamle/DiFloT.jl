"""
    solve_quickest_flow(net, s, t, F; method="INT-SG", solver="HiGHS", solver_method="primal", max_iter=100, tol=1e-13)

Solves the Quickest Flow Problem (QFP).

# Arguments
- `net`: The `NetworkInstance`.
- `s`, `t`: Source and sink nodes.
- `F`: The demand (flow value).
- `method`: `"BS"`, `"IBS"`, `"INT"`, `"INT-SG"`, `"NEWT"`, or `"LP"`.
- `solver`: `"HiGHS"` or `"Gurobi"`.
- `solver_method`: `"primal"`, `"dual"`, or `"network"`.
- `max_iter`, `tol`: Parameters for the inner root-finding in interpolation methods.
"""
function solve_quickest_flow(net, s::Int, t::Int, F::Number; method::String="INT-SG", solver::String="HiGHS", solver_method::String="primal", max_iter::Int=100, tol::Float64=1e-7)
    if method == "BS"
        return solve_qf_bs(net, s, t, F; solver=solver, solver_method=solver_method)
    elseif method == "IBS"
        return solve_qf_ibs(net, s, t, F; solver=solver, solver_method=solver_method)
    elseif method == "INT"
        return solve_qf_intp(net, s, t, F; solver=solver, solver_method=solver_method, max_iter=max_iter, tol=tol)
    elseif method == "INT-SG"
        return solve_qf_intsg(net, s, t, F; solver=solver, solver_method=solver_method, max_iter=max_iter, tol=tol)
    elseif method == "NEWT"
        return solve_qf_newton(net, s, t, F; solver=solver, solver_method=solver_method)
    elseif method == "LP"
        return solve_qf_lp(net, s, t, F; solver=solver, solver_method=solver_method)
    else
        error("Unknown method '$method'.")
    end
end
