# Interpolation search for the QFP (pure and safeguarded variants).

"""
    _illinois(P, p_lo, p_hi, max_iter, tol; x_tol=tol)

Bracketed Illinois method on [0,1]. Only a repeatedly retained endpoint
has its stored residual halved. The returned root is approximate.
"""
function _illinois(P, p_lo::Float64, p_hi::Float64, max_iter::Int, tol::Float64; x_tol::Float64=tol)
    lo, hi = 0.0, 1.0
    p_lo == 0.0 && return lo
    p_hi == 0.0 && return hi
    p_lo < 0.0 < p_hi || return NaN
    
    # Track which endpoint was replaced in the preceding iteration.
    last = 0
    for _ in 1:max_iter
        hi - lo <= x_tol && return (lo + hi) / 2

        # Regula falsi step; fall back to bisection if necessary.
        x = lo - p_lo * (hi - lo) / (p_hi - p_lo)
        if !(lo < x < hi)
            x = (lo + hi) / 2
        end
    
        px = P(x)
        isfinite(px) || return NaN
        abs(px) <= tol && return x

        # Halve the residual of an endpoint retained twice in succession.
        if px < 0.0
            lo, p_lo = x, px
            last == -1 && (p_hi *= 0.5)
            last = -1
        else
            hi, p_hi = x, px
            last == 1 && (p_lo *= 0.5)
            last = 1
        end
    end
    
    return (lo + hi) / 2
end

"""
    _interp_candidate(T_l,T_u,r_l,r_u,F,max_iter,tol)

Construct the interpolant of w(T)=v(T)-F. In normalized coordinates
q=(T_u-T)/L, its negative is
    P(q)=F-v_u+g_u*L*q-(g_u*L-(v_u-v_l))*q^b.
This avoids separately forming a and (T_u-T)^b. The projection guarantees
an interior integer test point even if the numerical root is inaccurate.
"""
function _interp_candidate(
    T_l::Int, T_u::Int, r_l::MCCPResult, r_u::MCCPResult, F::Int, max_iter::Int, tol::Float64
)
    L = T_u - T_l
    v_l, v_u = _value(r_l), _value(r_u)
    g_l, g_u = _slope(r_l), _slope(r_u)
    rise = v_u - v_l
    
    # Exact integer arithmetic for the coefficient and the degenerate case.
    D = Int128(g_u) * L - rise
    B = Int128(g_u - g_l) * L
    
    if D < 0 || B < D || g_l < 1 || !(v_l < F <= v_u)
        error("Inconsistent interpolation bracket or subgradients.")
    end
    
    if D == 0 || B == D
        # Linear interpolant: compute its ceiling exactly.
        r_ceil = T_l + Int(cld(Int128(L) * (F - v_l), Int128(rise)))
        return clamp(r_ceil, T_l + 1, T_u - 1)
    end

    # Nonlinear interpolant in normalized coordinates.
    b = Float64(B) / Float64(D)
    c = Float64(D)
    guL = Float64(Int128(g_u) * L)
    z_u = Float64(F - v_u)
    P(q) = muladd(guL, q, z_u) - c * q^b
    
    # Use endpoint values directly to avoid cancellation at q=1.
    # Root-finder width tolerance is in units of time-horizon.
    q = _illinois(P, z_u, Float64(F - v_l), max_iter, tol; x_tol=tol / Float64(L))
    isfinite(q) || return _midpoint(T_l, T_u)
    
    # ceil(T_u-L*q) = T_u-floor(L*q). Avoid cancellation of large horizons.
    T_cand = T_u - floor(Int, clamp(q, 0.0, 1.0) * L)
    return clamp(T_cand, T_l + 1, T_u - 1)
end

"""Shared implementation of pure and safeguarded interpolation search."""
function _solve_qf_interpolation(
    net, s::Int, t::Int, F::Real, safeguarded::Bool;
    solver::String="HiGHS", solver_method::String="primal", max_flow_solver::String="GraphsFlows",
    max_iter::Int=100, tol::Float64=1e-7, with_paths::Bool=false, verify::Bool=false
)
    F_int = _positive_demand(F)
    max_iter > 0 && isfinite(tol) && tol > 0 || throw(ArgumentError("Invalid root-finding parameters."))
    
    # Initialize the lower bound and its MCCP solution.
    init = _initial_lower(
        net, s, t, F_int;
        solver=solver, solver_method=solver_method, max_flow_solver=max_flow_solver
    )
    isnothing(init) && return _infeasible_qf()
    (; session, T_l, r_l) = init

    # The lower bound is optimal if it already satisfies the demand.
    if _value(r_l) >= F_int
        return _qf_result(
            session, T_l, r_l, F_int;
            with_paths=with_paths, verify=verify
        )
    end

    # Establish the bracket v(T_l) < F <= v(T_u).
    T_u, r_u = _initial_upper(session, T_l, r_l, F_int)
    activations = 0

    while T_l + 1 < T_u
        # Interpolation candidate (rounded and projected) and tentative bracket length.
        L = T_u - T_l
        T_cand = _interp_candidate(T_l, T_u, r_l, r_u, F_int, max_iter, tol)
        r_cand = solve_mccp(session, T_cand)
        
        feasible = _value(r_cand) >= F_int
        
        L_new = feasible ? T_cand - T_l : T_u - T_cand

        # Safeguard: require the interval to shrink by at least 1/5.
        if safeguarded && 5 * Int128(L_new) > 4 * Int128(L)
            # Midpoint of the original interval, not the tentative interval.
            activations += 1
            T_m = _midpoint(T_l, T_u)
            r_m = T_m == T_cand ? r_cand : solve_mccp(session, T_m)
        
            # Retain both oracle evaluations, exactly as in the paper.
            if feasible
                T_u, r_u = T_cand, r_cand
            else
                T_l, r_l = T_cand, r_cand
            end

            # Tighten the updated bracket using the midpoint solution.
            if _value(r_m) >= F_int
                if T_m < T_u
                    T_u, r_u = T_m, r_m
                end
            else
                if T_m > T_l
                    T_l, r_l = T_m, r_m
                end
            end
        elseif feasible
            T_u, r_u = T_cand, r_cand
        else
            T_l, r_l = T_cand, r_cand
        end
        
        T_l < T_u && _value(r_l) < F_int <= _value(r_u) || error("Invalid interpolation bracket.")
    end
    return _qf_result(
        session, T_u, r_u, F_int;
        with_paths=with_paths, verify=verify, predecessor_checked=true, safeguard_activations=activations
    )
end

"""Quickest Flow via pure interpolation search."""
solve_qf_intp(net, s::Int, t::Int, F::Real; kwargs...) =
    _solve_qf_interpolation(net, s, t, F, false; kwargs...)

"""Quickest Flow via safeguarded interpolation search."""
solve_qf_intsg(net, s::Int, t::Int, F::Real; kwargs...) =
    _solve_qf_interpolation(net, s, t, F, true; kwargs...)
