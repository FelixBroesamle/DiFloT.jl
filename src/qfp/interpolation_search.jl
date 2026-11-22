# =============================================================================
# Interpolation search for the Quickest Flow Problem
# =============================================================================
# Implements Algorithms 6 (pure) and 7 (safeguarded).
# =============================================================================

"""
    _eval_float(session, T) -> (time, v, g, res)

Evaluate the MCCP oracle at an integer horizon `T`. Returns a named tuple
with float fields for the interpolant arithmetic and the raw `MCCPResult`.
"""
function _eval_float(session::MCCPSession, T::Int)
    r = solve_mccp(session, T)
    v = r.status == STATUS_OPTIMAL ? Float64(r.v) : 0.0
    g = r.status == STATUS_OPTIMAL ? Float64(r.g) : 0.0
    return (time = Float64(T), v = v, g = g, res = r)
end

"""
    _interp_root(a, b, g_u, z_u, xmax, max_iter, tol) -> Float64

Find the unique root of `P(x) = g_u * x + a * x^b + z_u` on `(0, xmax)`
using the Illinois method.

`P` is convex and `P(0) = z_u < 0`, `P(xmax) ≥ 0` (guaranteed by the
sandwich lemma), so a bracketed root exists. Falls back to the tangent
root `-z_u/g_u` when the bracket is degenerate.
"""
function _interp_root(a::Float64, b::Float64, g_u::Float64, z_u::Float64, xmax::Float64, max_iter::Int, tol::Float64)
    P(x) = g_u * x + a * x^b + z_u
    x_lo = 1e-12
    x_hi = xmax
    P_lo = P(x_lo)
    P_hi = P(x_hi)
    
    P_lo == 0.0 && return x_lo
    P_hi == 0.0 && return x_hi
    P_lo * P_hi > 0 && return abs(P_lo) < abs(P_hi) ? x_lo : x_hi

    local x_new, P_new
    for _ in 1:max_iter
        (x_hi - x_lo) < tol && break
        x_new = x_hi - P_hi * (x_hi - x_lo) / (P_hi - P_lo)
        if !(x_lo < x_new < x_hi)
            x_new = 0.5 * (x_lo + x_hi)
        end
        P_new = P(x_new)
        abs(P_new) < tol && return x_new
        if P_new * P_hi < 0
            x_lo = x_new; P_lo = P_new
            P_hi *= 0.5
        else
            x_hi = x_new; P_hi = P_new
            P_lo *= 0.5
        end
    end
    return 0.5 * (x_lo + x_hi)
end

"""
    _interp_candidate(session, T_l, T_u, eval_l, eval_u, F_float, max_iter, tol)
        -> (T_cand::Int, eval_cand)

Compute the interpolant root on `[T_l, T_u]` and return the ceiled,
clamped candidate `T_cand ∈ {T_l+1, …, T_u-1}` together with its oracle
evaluation.
"""
function _interp_candidate(session::MCCPSession, T_l::Int64, T_u::Int64,
                           eval_l, eval_u, F_float::Float64,
                           max_iter::Int, tol::Float64)
    L = T_u - T_l
    v_l, g_l = eval_l.v, eval_l.g
    v_u, g_u = eval_u.v, eval_u.g
    z_l = F_float - v_l
    z_u = F_float - v_u
    secant_slope = (z_u - z_l) / L

    # ---------------------------------------------------------------------
    # Step 1 ("Compute the root r ∈ [T_l, T_u] of ŵ(T) = 0")
    # ---------------------------------------------------------------------
    # The interpolant is parameterised in shifted coordinates x = T_u − T:
    #     ŵ(T) = P(T_u − T) with P(x) = g_u·x + a·x^b + z_u.
    # So the root in T-space is r = T_u − x_root, where P(x_root) = 0.
    x_root = if abs(g_u + secant_slope) < 1e-9
        -z_u / g_u                                   # linear (tangent) case
    else
        b = (g_u - g_l) / (g_u + secant_slope)
        if abs(b) < 1e-9 || !isfinite(b)
            -z_u / g_u
        else
            a = -(g_u + secant_slope) / L^(b - 1)
            _interp_root(a, b, g_u, z_u, Float64(L), max_iter, tol)
        end
    end
    r = T_u - x_root                                 

    # ---------------------------------------------------------------------
    # Step 2 ("T_m ← max{T_l + 1, min{T_u − 1, ⌈r⌉}}")
    # ---------------------------------------------------------------------
    T_cand    = max(T_l + 1, min(T_u - 1, ceil(Int64, r)))
    eval_cand = _eval_float(session, Int(T_cand))
    return T_cand, eval_cand
end

"""
    solve_qf_intp(net, s, t, F; solver="HiGHS", solver_method="primal", max_iter=100, tol=1e-9)

Quickest Flow via Pure Interpolation Search.
"""
function solve_qf_intp(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal", max_iter::Int=100, tol::Float64=1e-9)
    F_int   = round(Int64, F)
    F_float = Float64(F_int)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    T0 = round(Int64, T0f)

    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) &&
        return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = round(Int64, mf.value)

    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)

    # Initial bracket
    eval0 = _eval_float(session, 0)
    T_l = max(T0, max(0, cld(F_int - round(Int64, eval0.v), gmax)))
    eval_l = _eval_float(session, Int(T_l))
    round(Int64, eval_l.v) >= F_int &&
        return QFResult(STATUS_OPTIMAL, Float64(T_l), eval_l.res.flows)

    g_l0 = let g = round(Int64, eval_l.g); g > 0 ? g : gmax end
    T_u = T_l + cld(F_int - round(Int64, eval_l.v), g_l0)
    eval_u = _eval_float(session, Int(T_u))

    # Main loop
    while T_l + 1 < T_u
        T_m, eval_m = _interp_candidate(session, T_l, T_u, eval_l, eval_u,
                                        F_float, max_iter, tol)
        if round(Int64, eval_m.v) >= F_int
            T_u, eval_u = T_m, eval_m
        else
            T_l, eval_l = T_m, eval_m
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(T_u), eval_u.res.flows)
end

"""
    solve_qf_intsg(net, s, t, F; solver="HiGHS", solver_method="primal", max_iter=100, tol=1e-9)

Quickest Flow via Safeguarded Interpolation Search.
"""
function solve_qf_intsg(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal", max_iter::Int=100, tol::Float64=1e-9)
    F_int   = round(Int64, F)
    F_float = Float64(F_int)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    T0 = round(Int64, T0f)

    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) &&
        return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = round(Int64, mf.value)

    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)

    # Initial bracket (Lemma 2.9)
    eval0 = _eval_float(session, 0)
    T_l = max(T0, max(0, cld(F_int - round(Int64, eval0.v), gmax)))
    eval_l = _eval_float(session, Int(T_l))
    round(Int64, eval_l.v) >= F_int &&
        return QFResult(STATUS_OPTIMAL, Float64(T_l), eval_l.res.flows)

    g_l0 = let g = round(Int64, eval_l.g); g > 0 ? g : gmax end
    T_u = T_l + cld(F_int - round(Int64, eval_l.v), g_l0)
    eval_u = _eval_float(session, Int(T_u))

    # Main loop — Algorithm 7
    while T_l + 1 < T_u
        L = T_u - T_l

        # --- Evaluate interpolant candidate ---
        T_cand, eval_cand = _interp_candidate(session, T_l, T_u, eval_l, eval_u,
                                              F_float, max_iter, tol)

        # --- Compute the length that WOULD result from accepting T_cand ---
        L_new = round(Int64, eval_cand.v) >= F_int ? (T_cand - T_l) : (T_u - T_cand)

        if 5 * L_new <= 4 * L
            # ---- Accept candidate ----
            if round(Int64, eval_cand.v) >= F_int
                T_u, eval_u = T_cand, eval_cand
            else
                T_l, eval_l = T_cand, eval_cand
            end
        else
            # ---- Safeguard: fall back to midpoint ----
            T_m = T_l + div(L, 2)
            eval_m = _eval_float(session, Int(T_m))
            if round(Int64, eval_m.v) >= F_int
                T_u, eval_u = T_m, eval_m
            else
                T_l, eval_l = T_m, eval_m
            end
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(T_u), eval_u.res.flows)
end
