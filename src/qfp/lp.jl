"""
    solve_qf_lp(net, s, t, F; solver="HiGHS", solver_method="primal")

Quickest Flow via the Lin-Jaillet parametric LP formulation.
"""
function solve_qf_lp(net, s::Int, t::Int, F::Number; solver::String="HiGHS", solver_method::String="primal")
    model = Model(_optimizer(solver; solver_method=solver_method, warm_start=false))
    set_silent(model)

    x = Dict{Tuple{Int,Int,Int},JuMP.VariableRef}()
    for e in net.edges
        k = arc_key(e)
        x[k] = @variable(model, lower_bound=0.0)
    end
    theta = @variable(model, lower_bound=0.0)

    # The objective is exactly the continuous time horizon T_cont*
    @objective(model, Min, F * theta + sum(e.transit_time * x[arc_key(e)] for e in net.edges))

    in_keys, out_keys = incidence_mapping(net)
    for v in 1:net.num_nodes
        rhs = v == s ? 1.0 : (v == t ? -1.0 : 0.0)
        @constraint(model, sum(x[k] for k in out_keys[v]; init=0.0) - sum(x[k] for k in in_keys[v]; init=0.0) == rhs)
    end

    for e in net.edges
        @constraint(model, x[arc_key(e)] <= e.up * theta)
    end

    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL
        return QFResult(STATUS_INFEASIBLE, -1.0, Float64[])
    end

    # THE FIX: The optimal objective value IS the continuous time horizon tau*
    tau_star = objective_value(model)
    
    # Round to the discrete time horizon
    T = max(0, ceil(Int, tau_star) - 1)
    
    session = MCCPSession(net, s, t; solver=solver, solver_method=solver_method)
    r = solve_mccp(session, T)
    
    # Oracle verification to handle floating point noise
    if r.v < F
        r = solve_mccp(session, T + 1)
        return QFResult(STATUS_OPTIMAL, Float64(T + 1), r.flows)
    elseif T >= 1 && solve_mccp(session, T - 1).v >= F
        r_prev = solve_mccp(session, T - 1)
        return QFResult(STATUS_OPTIMAL, Float64(T - 1), r_prev.flows)
    else
        return QFResult(STATUS_OPTIMAL, Float64(T), r.flows)
    end
end