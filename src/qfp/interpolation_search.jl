function _eval_float(session::MCCPSession, T::Int)
    r = solve_mccp(session, T)
    v = r.status == STATUS_OPTIMAL ? Float64(r.v) : 0.0
    g = r.status == STATUS_OPTIMAL ? Float64(r.g) : 0.0
    return (time = Float64(T), v = v, g = g, res = r)
end

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
    solve_qf_intp(net, s, t, F; solver="HiGHS", solver_method="primal", max_iter=100, tol=1e-13)

Quickest Flow via Pure Interpolation Search.
"""
function solve_qf_intp(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal", max_iter::Int=100, tol::Float64=1e-13)
    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)
    F_target = Float64(F)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    T0 = round(Int, T0f)

    mf = ORToolsBackend.solve_max_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = Float64(mf.value)

    eval0 = _eval_float(session, 0)
    T_l = max(T0, ceil(Int, max(0.0, (F_target - eval0.v) / gmax)))
    eval_l = _eval_float(session, T_l)
    eval_l.v >= F_target && return QFResult(STATUS_OPTIMAL, Float64(eval_l.time), eval_l.res.flows)

    g_l = eval_l.g > 0 ? eval_l.g : gmax
    T_u = T_l + ceil(Int, (F_target - eval_l.v) / g_l)
    eval_u = _eval_float(session, T_u)

    iter = 0
    while T_l + 1 < T_u
        iter += 1
        if iter % 10 == 0
            @info "INT: iter=$iter, T_l=$T_l, T_u=$(T_u)"
        end

        v_l, g_l = eval_l.v, eval_l.g
        v_u, g_u = eval_u.v, eval_u.g
        z_l, z_u = F_target - v_l, F_target - v_u
        secant_slope = (z_u - z_l) / (T_u - T_l)

        x_root = 0.0
        if abs(g_u + secant_slope) < 1e-9
            x_root = -z_u / g_u
        else
            b = (g_u - g_l) / (g_u + secant_slope)
            if abs(b) < 1e-9 || !isfinite(b)
                x_root = -z_u / g_u
            else
                a = -(g_u + secant_slope) / (T_u - T_l)^(b - 1)
                x_root = _interp_root(a, b, g_u, z_u, Float64(T_u - T_l), max_iter, tol)
            end
        end

        T_m = clamp(ceil(Int, T_u - x_root), T_l + 1, T_u - 1)
        eval_m = _eval_float(session, T_m)
        
        if eval_m.v >= F_target
            T_u, eval_u = T_m, eval_m
        else
            T_l, eval_l = T_m, eval_m
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(T_u), eval_u.res.flows)
end

"""
    solve_qf_intsg(net, s, t, F; solver="HiGHS", solver_method="primal", max_iter=100, tol=1e-13)

Quickest Flow via Safeguarded Interpolation Search.
"""
function solve_qf_intsg(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal", max_iter::Int=100, tol::Float64=1e-13)
    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)
    F_target = Float64(F)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    T0 = round(Int, T0f)

    mf = ORToolsBackend.solve_max_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = Float64(mf.value)

    eval0 = _eval_float(session, 0)
    T_l = max(T0, ceil(Int, max(0.0, (F_target - eval0.v) / gmax)))
    eval_l = _eval_float(session, T_l)
    eval_l.v >= F_target && return QFResult(STATUS_OPTIMAL, Float64(eval_l.time), eval_l.res.flows)

    g_l = eval_l.g > 0 ? eval_l.g : gmax
    T_u = T_l + ceil(Int, (F_target - eval_l.v) / g_l)
    eval_u = _eval_float(session, T_u)

    iter = 0
    while T_l + 1 < T_u
        iter += 1
        if iter % 10 == 0
            @info "INT-SG: iter=$iter, T_l=$T_l, T_u=$(T_u)"
        end

        v_l, g_l = eval_l.v, eval_l.g
        v_u, g_u = eval_u.v, eval_u.g
        z_l, z_u = F_target - v_l, F_target - v_u
        secant_slope = (z_u - z_l) / (T_u - T_l)

        x_root = 0.0
        if abs(g_u + secant_slope) < 1e-9
            x_root = -z_u / g_u
        else
            b = (g_u - g_l) / (g_u + secant_slope)
            if abs(b) < 1e-9 || !isfinite(b)
                x_root = -z_u / g_u
            else
                a = -(g_u + secant_slope) / (T_u - T_l)^(b - 1)
                x_root = _interp_root(a, b, g_u, z_u, Float64(T_u - T_l), max_iter, tol)
            end
        end

        L = T_u - T_l
        T_cand = ceil(Int, T_u - x_root)
        T_min = T_l + ceil(Int, L / 4)
        T_max = T_u - ceil(Int, L / 4)
        T_m = max(T_min, min(T_max, T_cand))

        eval_m = _eval_float(session, T_m)
        
        if eval_m.v >= F_target
            T_u, eval_u = T_m, eval_m
        else
            T_l, eval_l = T_m, eval_m
        end
    end
    return QFResult(STATUS_OPTIMAL, Float64(T_u), eval_u.res.flows)
end