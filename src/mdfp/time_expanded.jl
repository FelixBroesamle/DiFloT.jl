function solve_max_dynamic_flow_te(net, s::Int, t::Int, T::Int)
    n = net.num_nodes
    starts, ends, caps = Int[], Int[], Int[]
    s_super = (T + 1) * n
    t_super = (T + 1) * n + 1
    c_rate = c_rate_bound(net, s, t)
    c_total = (T + 1) * c_rate

    for v in 1:n, theta in 0:(T-1)
        push!(starts, theta * n + (v - 1))
        push!(ends,   (theta + 1) * n + (v - 1))
        push!(caps,   c_total)
    end

    for e in net.edges
        for theta in 0:(T - e.transit_time)
            push!(starts, theta * n + (e.tail - 1))
            push!(ends,   (theta + e.transit_time) * n + (e.head - 1))
            push!(caps,   e.up)
        end
    end

    for theta in 0:T
        push!(starts, s_super); push!(ends, theta * n + (s - 1)); push!(caps, c_rate)
        push!(starts, theta * n + (t - 1)); push!(ends, t_super); push!(caps, c_rate)
    end

    smf = ORToolsBackend.pygraph_mf.SimpleMaxFlow()
    smf.add_arcs_with_capacity(ORToolsBackend.pynp.array(starts), ORToolsBackend.pynp.array(ends), ORToolsBackend.pynp.array(caps))
    status = smf.solve(s_super, t_super)
    
    if pyconvert(Bool, status == smf.OPTIMAL)
        val = pyconvert(Int64, smf.optimal_flow())
        return MDFResult(STATUS_OPTIMAL, Float64(val), Float64[], Tuple{Vector{Int},Float64}[])
    else
        return MDFResult(STATUS_INFEASIBLE, 0.0, Float64[], Tuple{Vector{Int},Float64}[])
    end
end

function solve_max_dynamic_flow_te_lp(net, s::Int, t::Int, T::Int, solver::String, solver_method::String)
    model = Model(_optimizer(solver; solver_method=solver_method))
    set_silent(model)
    
    n = net.num_nodes
    s_super = (T + 1) * n + 1
    t_super = (T + 1) * n + 2
    c_rate = c_rate_bound(net, s, t)
    c_total = (T + 1) * c_rate
    
    x = Dict{Tuple{Int,Int,Int}, JuMP.VariableRef}()
    
    for v in 1:n, theta in 0:(T-1)
        k = (theta * n + v, (theta + 1) * n + v, 0)
        x[k] = @variable(model, lower_bound=0.0, upper_bound=Float64(c_total))
    end
    
    for (i, e) in enumerate(net.edges)
        for theta in 0:(T - e.transit_time)
            k = (theta * n + e.tail, (theta + e.transit_time) * n + e.head, i)
            x[k] = @variable(model, lower_bound=0.0, upper_bound=Float64(e.up))
        end
    end
    
    for theta in 0:T
        x[(s_super, theta * n + s, 0)] = @variable(model, lower_bound=0.0, upper_bound=Float64(c_rate))
        x[(theta * n + t, t_super, 0)] = @variable(model, lower_bound=0.0, upper_bound=Float64(c_rate))
    end
    
    in_edges = Dict{Int, Vector{Tuple{Int,Int,Int}}}()
    out_edges = Dict{Int, Vector{Tuple{Int,Int,Int}}}()
    for k in keys(x)
        push!(get!(out_edges, k[1], Tuple{Int,Int,Int}[]), k)
        push!(get!(in_edges, k[2], Tuple{Int,Int,Int}[]), k)
    end
    
    for v in 1:(T+1)*n
        inflow = sum(x[k] for k in get(in_edges, v, Tuple{Int,Int,Int}[]); init=0.0)
        outflow = sum(x[k] for k in get(out_edges, v, Tuple{Int,Int,Int}[]); init=0.0)
        @constraint(model, inflow == outflow)
    end
    
    @objective(model, Max, sum(x[k] for k in get(out_edges, s_super, Tuple{Int,Int,Int}[])))
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        return MDFResult(STATUS_OPTIMAL, objective_value(model), Float64[], Tuple{Vector{Int},Float64}[])
    else
        return MDFResult(STATUS_INFEASIBLE, 0.0, Float64[], Tuple{Vector{Int},Float64}[])
    end
end