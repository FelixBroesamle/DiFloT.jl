"""
    solve_qf_bs(net, s, t, F; solver="HiGHS", solver_method="primal")

Quickest Flow via standard binary search.
"""
function solve_qf_bs(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    T0 = Int(T0f)

    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = Float64(mf.value)

    r0 = solve_mccp(session, 0)
    v0 = r0.status == STATUS_OPTIMAL ? r0.v : 0.0

    Tl = max(T0, ceil(Int, max(0.0, (F - v0) / gmax)))
    rl = solve_mccp(session, Tl)
    rl.v >= F && return QFResult(STATUS_OPTIMAL, Float64(Tl), rl.flows)

    gl = rl.g > 0 ? rl.g : gmax
    Tu = Tl + ceil(Int, (F - rl.v) / gl)
    ru = solve_mccp(session, Tu)

    iter = 0
    while Tl + 1 < Tu
        iter += 1
        if iter % 10 == 0
            @info "BS: iter=$iter, T_l=$Tl, T_u=$Tu"
        end
        mid = div(Tl + Tu, 2)
        rm = solve_mccp(session, mid)
        if rm.v >= F
            Tu = mid
            ru = rm
        else
            Tl = mid
            rl = rm
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(Tu), ru.flows)
end
