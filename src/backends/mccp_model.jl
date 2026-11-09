mutable struct MCCPSession
    net::NetworkInstance
    s::Int
    t::Int
    model::JuMP.Model
    x::Dict{Tuple{Int,Int,Int},JuMP.VariableRef}
    x_ret::JuMP.VariableRef
    last_coefficient::Float64
    last_result::Union{Nothing,MCCPResult}
end

"""
    MCCPSession(net, s, t; solver="HiGHS", solver_method="primal")

Builds the parametric circulation model once for warm-starting. 
"""
function MCCPSession(net::NetworkInstance, s::Int, t::Int; solver::String="HiGHS", solver_method::String="primal")
    model = Model(_optimizer(solver; solver_method=solver_method, warm_start=true))
    set_silent(model)

    x = Dict{Tuple{Int,Int,Int},JuMP.VariableRef}()
    for e in net.arcs
        k = arc_key(e)
        x[k] = @variable(model, lower_bound=Float64(e.low), upper_bound=Float64(e.up))
    end

    c_rate = c_rate_bound(net, s, t)
    x_ret = @variable(model, lower_bound=0.0, upper_bound=Float64(c_rate))

    in_keys, out_keys = incidence_mapping(net)
    for v in 1:net.num_nodes
        inflow = sum(x[k] for k in in_keys[v]; init=0.0)
        outflow = sum(x[k] for k in out_keys[v]; init=0.0)
        if v == s; inflow += x_ret; end
        if v == t; outflow += x_ret; end
        @constraint(model, inflow == outflow)
    end

    fixed_obj = sum(e.transit_time * x[arc_key(e)] for e in net.arcs)
    @objective(model, Min, fixed_obj + 0.0 * x_ret)

    return MCCPSession(net, s, t, model, x, x_ret, NaN, nothing)
end

"""
    solve_mccp(session::MCCPSession, T::Int; tol=1e-9, force=false) -> MCCPResult

Solves the parametric minimum cost circulation for a given time horizon `T`.
"""
function solve_mccp(session::MCCPSession, T::Int; tol::Real=1e-9, force::Bool=false)
    a = Float64(-(T + 1))
    if !force && session.last_coefficient == a && session.last_result !== nothing
        return session.last_result
    end

    JuMP.set_objective_coefficient(session.model, session.x_ret, a)
    session.last_coefficient = a
    JuMP.optimize!(session.model)

    st = termination_status(session.model)
    if st != MOI.OPTIMAL
        session.last_result = MCCPResult(STATUS_INFEASIBLE, -Inf, 0.0, Float64[])
        return session.last_result
    end

    flows = Vector{Float64}(undef, session.net.num_arcs)
    for (i, e) in enumerate(session.net.arcs)
        v = JuMP.value(session.x[arc_key(e)])
        flows[i] = abs(v) < tol ? 0.0 : v
    end

    session.last_result = MCCPResult(STATUS_OPTIMAL,
        -objective_value(session.model),
        JuMP.value(session.x_ret),
        flows)
    return session.last_result
end
