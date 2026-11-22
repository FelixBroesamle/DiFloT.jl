# =============================================================================
# Newton-type search for the Quickest Flow Problem
# =============================================================================
# Implements Algorithm 5 from the paper (§3.3), with:
#   - Integer-demand enforcement
#   - Exact T_l via rounded g_max
#   - Safety cap on iterations
#   - Floor-then-+1 verification (matches paper's two-step check)
# =============================================================================

"""
    _eval_param(session::MCCPSession, theta::Float64) -> (v, g)

Evaluate the MCCP oracle at a real-valued `theta`. Returns `(v, g)`.
On solver failure, returns `(0.0, 0.0)` — the caller must guard.
"""
function _eval_param(session::MCCPSession, theta::Float64)
    coeff = -(theta + 1.0)
    JuMP.set_objective_coefficient(session.model, session.x_ret, coeff)
    session.last_coefficient = coeff
    JuMP.optimize!(session.model)

    if termination_status(session.model) != MOI.OPTIMAL
        return (v = 0.0, g = 0.0)
    end

    flows = Vector{Float64}(undef, session.net.num_arcs)
    for (i, e) in enumerate(session.net.arcs)
        xv = JuMP.value(session.x[arc_key(e)])
        flows[i] = abs(xv) < 1e-9 ? 0.0 : xv
    end
    session.last_result = MCCPResult(
        STATUS_OPTIMAL,
        -objective_value(session.model),
        JuMP.value(session.x_ret),
        flows,
    )
    return (v = session.last_result.v, g = session.last_result.g)
end

"""
    solve_qf_newton(net, s, t, F; solver="HiGHS", solver_method="primal")

Solve the Quickest Flow Problem via the Newton-type search of Burkard et al.
This implements Algorithm 5.

The demand `F` must be integer-valued; otherwise the termination criterion
`|v(θ) - F| < 1` no longer implies `T* ∈ {⌊θ⌋, ⌊θ⌋+1}`.
"""
function solve_qf_newton(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
    Fq = Float64(F)

    # ----------------------------------------------------------------- 
    # Certificates
    # -----------------------------------------------------------------
    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    
    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    
    gmax = round(Int64, mf.value)

    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)

    # -----------------------------------------------------------------
    # Initial lower bound
    # -----------------------------------------------------------------
    e0 = _eval_param(session, 0.0)
    Tl = max(round(Int64, T0f), ceil(Int64, max(0.0, (Fq - e0.v) / gmax)))
    theta = Float64(Tl)
    cur = _eval_param(session, theta)

    # Early return: T_l already feasible
    if cur.v >= Fq
        return QFResult(STATUS_OPTIMAL, Float64(Tl), session.last_result.flows)
    end

    # -----------------------------------------------------------------
    # Newton iterations: paper's criterion `While |v(θ) - F| ≥ 1`
    # -----------------------------------------------------------------
    iter = 0
    max_iter = 1000
    while abs(cur.v - Fq) >= 1.0
        iter += 1
        iter > max_iter && error(
            "Newton failed to converge in $max_iter iterations " *
            "(θ = $theta, v = $(cur.v), g = $(cur.g), F = $Fq)"
        )

        # Safety: slope must be positive to take a Newton step
        if cur.g <= 1e-9
            break
        end

        theta = theta + (Fq - cur.v) / cur.g
        cur   = _eval_param(session, theta)
    end

    # -----------------------------------------------------------------
    # Two-step verification: T* ∈ {⌊θ⌋, ⌊θ⌋+1}
    # -----------------------------------------------------------------
    T = floor(Int64, theta)
    r = _eval_param(session, Float64(T))
    if r.v >= Fq
        return QFResult(STATUS_OPTIMAL, Float64(T), session.last_result.flows)
    end

    _eval_param(session, Float64(T + 1))
    return QFResult(STATUS_OPTIMAL, Float64(T + 1), session.last_result.flows)
end
