# Shared construction for combinatorial and LP-based time expansion.
function _time_expanded_arcs(net, s::Int, t::Int, T::Int)
    _validate_network(net, s, t)
    _check_horizon(net, s, t, T)

    # Node v at time theta has index theta*n + v.
    n = net.num_nodes
    N = Base.Checked.checked_add(Base.Checked.checked_mul(T + 1, n), 2)
    s_super, t_super = N - 1, N
    c_rate = c_rate_bound(net, s, t)
    c_total = (T + 1) * c_rate

    # Arc arrays: tail, head, and capacity.
    I, J, V = Int[], Int[], Int[]
    function add_arc!(u, v, cap)
        cap > 0 || return
        push!(I, u)
        push!(J, v)
        push!(V, cap)
        return nothing
    end

    # Holdover arcs; c_total is a valid finite replacement for infinity.
    for theta in 0:T-1, v in 1:n
        add_arc!(theta * n + v, (theta + 1) * n + v, c_total)
    end

    # Movement arcs: include only departures arriving by T.
    for e in net.arcs
        # Self-loops cannot improve a maximum flow with holdover arcs.
        e.up > 0 && e.tail != e.head || continue
        for theta in 0:T-e.transit_time
            add_arc!(theta * n + e.tail, (theta + e.transit_time) * n + e.head, e.up)
        end
    end
    
    # Per-time terminal capacity c_rate preserves an optimal temporally
    # repeated flow and therefore the optimum of the expanded network.
    for theta in 0:T
        add_arc!(s_super, theta * n + s, c_rate)
        add_arc!(theta * n + t, t_super, c_rate)
    end
    
    return I, J, V, N, s_super, t_super, c_total
end

"""
    solve_max_dynamic_flow_te(net,s,t,T;solver="GraphsFlows")

Compute the MDFP value via a maximum flow on the time-expanded network.
"""
function solve_max_dynamic_flow_te(net, s::Int, t::Int, T::Int; solver::String="GraphsFlows")
    I, J, V, N, ss, tt, c_total = _time_expanded_arcs(net, s, t, T)
    f = isempty(V) ? 0 : _max_flow_value(I, J, V, N, ss, tt, solver)
    0 <= f <= c_total || error("Invalid time-expanded maximum-flow value.")
    
    return MDFResult(STATUS_OPTIMAL, Float64(f), Float64[], PathFlow[])
end

"""Internal LP version of time expansion, retained for compatibility."""
function solve_max_dynamic_flow_te_lp(net, s::Int, t::Int, T::Int, solver::String, solver_method::String)
    I, J, V, N, ss, tt, c_total = _time_expanded_arcs(net, s, t, T)
    model = Model(_optimizer(solver; solver_method=solver_method))
    set_silent(model)
    
    x = JuMP.VariableRef[]
    incoming = [Int[] for _ in 1:N]
    outgoing = [Int[] for _ in 1:N]
    
    for i in eachindex(V)
        push!(x, @variable(model, lower_bound=0.0, upper_bound=Float64(V[i])))
        push!(outgoing[I[i]], i)
        push!(incoming[J[i]], i)
    end
    
    for v in 1:N-2
        @constraint(
            model, 
            sum(x[i] for i in incoming[v]; init=0.0) == sum(x[i] for i in outgoing[v]; init=0.0)
        )
    end
    
    @objective(model, Max, sum(x[i] for i in outgoing[ss]; init=0.0))
    optimize!(model)
    _require_optimal(model, "Time-expanded LP")
    f = _rounded_flow(objective_value(model), "time-expanded objective")
    0 <= f <= c_total || error("Invalid time-expanded LP value.")
    
    return MDFResult(STATUS_OPTIMAL, Float64(f), Float64[], PathFlow[])
end
