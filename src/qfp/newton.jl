function _eval_param(session::MCCPSession, theta::Float64)
    JuMP.set_objective_coefficient(session.model, session.x_ret, -(theta + 1.0))
    session.last_coefficient = -(theta + 1.0)
    JuMP.optimize!(session.model)

    if termination_status(session.model) == MOI.OPTIMAL
        flows = Vector{Float64}(undef, session.net.num_arcs)
        for (i, e) in enumerate(session.net.arcs)
            v = JuMP.value(session.x[arc_key(e)])
            flows[i] = abs(v) < 1e-9 ? 0.0 : v
        end
        session.last_result = MCCPResult(STATUS_OPTIMAL, -objective_value(session.model), JuMP.value(session.x_ret), flows)
        return (v = session.last_result.v, g = session.last_result.g)
    else
        return (v = 0.0, g = 0.0)
    end
end

"""
    solve_qf_newton(net, s, t, F; solver="HiGHS", solver_method="primal")

Quickest Flow via Newton-type search.
"""
function solve_qf_newton(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
    Fq = Float64(F)

    T0f = compute_shortest_path(net, s, t)
    isinf(T0f) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    
    mf = compute_max_static_flow(net, s, t)
    (mf.status != STATUS_OPTIMAL || mf.value == 0) && return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    gmax = Float64(mf.value)

    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)

    e0 = _eval_param(session, 0.0)
    Tl = max(round(Int, T0f), ceil(Int, max(0.0, (Fq - e0.v) / gmax)))
    theta = Float64(Tl)
    cur = _eval_param(session, theta)

    cur.v >= Fq && return QFResult(STATUS_OPTIMAL, Float64(Tl), session.last_result.flows)

    iter = 0
    while true
        iter += 1
        if iter % 10 == 0
            @info "NEWT: iter=$iter, theta=$theta, v=$cur.v"
        end
        cur.g <= 1e-9 && return QFResult(STATUS_OPTIMAL, Float64(ceil(Int, theta)), session.last_result.flows)

        theta_next = theta + (Fq - cur.v) / cur.g

        if theta_next == theta || abs(cur.v - Fq) < 1.0
            T_up = ceil(Int, theta)
            prev = _eval_param(session, Float64(T_up - 1))
            if prev.v >= Fq
                return QFResult(STATUS_OPTIMAL, Float64(T_up - 1), session.last_result.flows)
            else
                _eval_param(session, Float64(T_up))
                return QFResult(STATUS_OPTIMAL, Float64(T_up), session.last_result.flows)
            end
        end

        theta = theta_next
        cur = _eval_param(session, theta)
    end
end
