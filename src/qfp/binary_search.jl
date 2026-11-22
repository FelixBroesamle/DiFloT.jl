"""
    solve_qf_bs(net, s, t, F; solver="HiGHS", solver_method="primal")

Quickest Flow via standard binary search.
"""
function solve_qf_bs(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
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

    while Tl + 1 < Tu
        mid = div(Tl + Tu, 2)
        rm = solve_mccp(session, mid)
        if round(Int64, rm.v) >= F_int
            Tu = mid
            ru = rm
        else
            Tl = mid
            rl = rm
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(Tu), ru.flows)
end
