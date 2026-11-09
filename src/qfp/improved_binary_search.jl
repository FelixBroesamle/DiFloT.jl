"""
    solve_qf_ibs(net, s, t, F; solver="HiGHS", solver_method="primal")

Quickest Flow via Improved Binary Search (chord + tangent bounds).
"""
function solve_qf_ibs(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)
    Fq = round(Int64, F)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    T0 = Int64(round(T0f))

    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = Int64(mf.value)

    r0 = solve_mccp(session, 0)
    Tl = max(T0, max(0, cld(Fq - round(Int64, r0.v), gmax)))
    rl = solve_mccp(session, Tl)
    rl.v >= Fq && return QFResult(STATUS_OPTIMAL, Float64(Tl), rl.flows)

    gl = rl.g > 0 ? round(Int64, rl.g) : gmax
    Tu = Tl + max(0, cld(Fq - round(Int64, rl.v), gl))
    ru = solve_mccp(session, Tu)

    iter = 0
    while Tl + 1 < Tu
        iter += 1
        if iter % 10 == 0
            @info "IBS: iter=$iter, T_l=$Tl, T_u=$Tu"
        end

        v_l, g_l = round(Int64, rl.v), round(Int64, rl.g)
        v_u, g_u = round(Int64, ru.v), round(Int64, ru.g)
        denom = Int128(v_u - v_l)
        alpha = denom > 0 ? Tl + Int64(cld(Int128(Tu - Tl) * Int128(Fq - v_l), denom)) : Tl + 1
        sigma_l = g_l > 0 ? Tl + cld(Fq - v_l, g_l) : typemax(Int64)
        sigma_u = g_u > 0 ? Tu + cld(Fq - v_u, g_u) : typemax(Int64)
        beta    = min(sigma_l, sigma_u)

        lo_c = max(Tl + 1, min(Tu - 1, alpha))
        hi_c = max(Tl + 1, min(Tu - 1, beta))
        Tm   = lo_c >= hi_c ? lo_c : div(lo_c + hi_c, 2)

        rm = solve_mccp(session, Tm)
        if rm.v >= Fq
            Tu, ru = Tm, rm
        else
            Tl, rl = Tm, rm
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(Tu), ru.flows)
end
