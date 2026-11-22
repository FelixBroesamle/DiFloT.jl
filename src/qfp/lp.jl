# =============================================================================
# Lin–Jaillet LP approach for the Quickest Flow Problem
# =============================================================================
# Implements Algorithm 8 of the paper, §3.5.
#
# The LP (Eq. 3.1) yields the continuous-time horizon τ*, whose rounding
# gives a candidate discrete horizon T within ±1 of the true optimum.
# Two exact MCCP oracle calls verify the candidate exactly, independent
# of the LP solver's numerical tolerance.
# =============================================================================

"""
    solve_qf_lp(net, s, t, F; solver="HiGHS", solver_method="primal")
        -> QFResult

Quickest Flow via the Lin–Jaillet parametric LP (Algorithm 8).

Solves one parametric LP (Eq. 3.1) and verifies the resulting candidate
horizon `T = max{0, ⌈τ̃⌉ − 1}` with at most two exact MCCP oracle calls.

# Arguments
- `net::NetworkInstance`: network. All arcs must satisfy `e.low = 0`
  (the classical Lin–Jaillet formulation); a non-zero lower bound
  raises `ArgumentError`.
- `s`, `t::Int`: source and sink.
- `F::Number`: demand; must be integer and `> 0`.
- `solver::String`: `"HiGHS"` (default) or `"Gurobi"`.
- `solver_method::String`: `"primal"` (default) or `"dual"`.

# Returns
- `QFResult(STATUS_OPTIMAL, T*, x_{T*})` on success.
- `QFResult(STATUS_INFEASIBLE, -1.0, [])` if `g_max = 0`.

# Complexity
1 LP solve + ≤ 2 MCCP oracle calls (Theorem 3.6, §3.5).
"""
function solve_qf_lp(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
    # ---- Input validation --------------------------------------------------
    isinteger(F) || throw(ArgumentError(
        "LP method requires integer demand F; got F = $F"))
    F > 0        || throw(ArgumentError(
        "LP method requires F > 0; got F = $F"))
    F_int = round(Int64, F)
    F_flt = Float64(F_int)

    any(e -> e.low > 0, net.arcs) &&
        throw(ArgumentError("LP method requires e.low = 0 on all arcs "))

    lowercase(solver) in ("highs", "gurobi") ||
        throw(ArgumentError(
            "LP method requires an LP-capable solver (\"HiGHS\" or \"Gurobi\"); " *
            "got solver=\"$solver\"."))

    # ---- feasibility certificate -------------------
    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) &&
        return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])

    
    model = Model(_optimizer(solver; solver_method=solver_method, warm_start=false))
    set_silent(model)

    x = Dict{Tuple{Int,Int,Int},JuMP.VariableRef}()
    for e in net.arcs
        k = arc_key(e)
        x[k] = @variable(model, lower_bound=0.0)
    end
    theta = @variable(model, lower_bound=0.0)

    # The objective is exactly the continuous time horizon T_cont*
    @objective(model, Min, F_flt * theta + sum(e.transit_time * x[arc_key(e)] for e in net.arcs))

    in_keys, out_keys = incidence_mapping(net)
    for v in 1:net.num_nodes
        rhs = v == s ? 1.0 : (v == t ? -1.0 : 0.0)
        @constraint(model, sum(x[k] for k in out_keys[v]; init=0.0) - sum(x[k] for k in in_keys[v]; init=0.0) == rhs)
    end

    for e in net.arcs
        @constraint(model, x[arc_key(e)] <= e.up * theta)
    end

    optimize!(model)
    termination_status(model) == MOI.OPTIMAL ||
        return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])

    tau_star = objective_value(model)
    (isfinite(tau_star) && tau_star >= 0) ||
        error("LP solver returned non-finite or negative objective value " *
              "tau* = $tau_star")
    
    T = max(0, ceil(Int64, tau_star) - 1)
    
    # ---- exact MCCP verification -----------------
    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)
    r = solve_mccp(session, Int(T))
    
    if round(Int64, r.v) < F_int
        r_next = solve_mccp(session, Int(T + 1))
        return QFResult(STATUS_OPTIMAL, Float64(T + 1), r_next.flows)
    end

    if T == 0
        return QFResult(STATUS_OPTIMAL, Float64(T), r.flows)
    else
        # single MCCP(T-1) call, stored in r_prev
        r_prev = solve_mccp(session, Int(T - 1))
        # If v_{T-1} ≥ F  
        if round(Int64, r_prev.v) >= F_int
            return QFResult(STATUS_OPTIMAL, Float64(T - 1), r_prev.flows)
        else
            return QFResult(STATUS_OPTIMAL, Float64(T), r.flows)
        end
    end
end
