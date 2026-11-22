"""
    solve_qf_ibs(net, s, t, F; solver="HiGHS", solver_method="primal")

Quickest Flow via Improved Binary Search (chord + tangent bounds).
"""
function solve_qf_ibs(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
    F_int = round(Int64, F)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    T0 = round(Int64, T0f)

    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) &&
        return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = Int64(mf.value)

    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)

    # Initial bracket
    r0 = solve_mccp(session, 0)
    v0_int = round(Int64, r0.v)
    Tl = max(T0, max(0, cld(F_int - v0_int, gmax)))
    rl = solve_mccp(session, Tl)
    rl_v_int = round(Int64, rl.v)
    rl_v_int >= F_int && return QFResult(STATUS_OPTIMAL, Float64(Tl), rl.flows)

    rl_g_int = round(Int64, rl.g)
    gl = rl_g_int > 0 ? rl_g_int : gmax
    Tu = Tl + cld(F_int - rl_v_int, gl)
    ru = solve_mccp(session, Tu)

    # Main loop
    while Tl + 1 < Tu
        v_l, g_l = round(Int64, rl.v), round(Int64, rl.g)
        v_u, g_u = round(Int64, ru.v), round(Int64, ru.g)

        # Secant lower bound (Int128 for overflow safety)
        denom = Int128(v_u - v_l)
        alpha = denom > 0 ? Tl + Int64(cld(Int128(Tu - Tl) * Int128(F_int - v_l), denom)) : Tl + 1

        # Tangent upper bounds
        sigma_l = g_l > 0 ? Tl + cld(F_int - v_l, g_l) : typemax(Int64)
        sigma_u = g_u > 0 ? Tu + cld(F_int - v_u, g_u) : typemax(Int64)
        beta    = min(sigma_l, sigma_u)

        lo_c = max(Tl + 1, min(Tu - 1, alpha))
        hi_c = max(Tl + 1, min(Tu - 1, beta))
        Tm   = lo_c >= hi_c ? lo_c : fld(lo_c + hi_c, 2) 

        rm = solve_mccp(session, Tm)
        if round(Int64, rm.v) >= F_int
            Tu, ru = Tm, rm
        else
            Tl, rl = Tm, rm
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(Tu), ru.flows)
end
