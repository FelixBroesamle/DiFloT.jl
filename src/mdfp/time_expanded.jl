"""
    solve_max_dynamic_flow_te(net, s, t, T)

Solves the MDFP using a Time-Expanded network and a highly optimized Sparse Matrix
with Dinic's Algorithm from `GraphsFlows.jl`.

# Arguments
- `net`: The `NetworkInstance`.
- `s`: The source node.
- `t`: The sink node.
- `T`: The time horizon.

# Returns
An `MDFResult`.
"""
function solve_max_dynamic_flow_te(net, s::Int, t::Int, T::Int)
    n = net.num_nodes
    s_super = (T + 1) * n + 1
    t_super = (T + 1) * n + 2
    num_nodes_te = t_super
    c_rate = c_rate_bound(net, s, t)
    c_total = Float64((T + 1) * c_rate)

    # Arrays for Sparse Matrix Construction
    I = Int[]
    J = Int[]
    V = Float64[]

    # Helper to add to sparse arrays
    function add_arc!(u, v, cap)
        push!(I, u); push!(J, v); push!(V, Float64(cap))
    end

    # 1. Holdover arcs
    for v in 1:n, theta in 0:(T-1)
        add_arc!(theta * n + v, (theta + 1) * n + v, c_total)
    end

    # 2. Network arcs
    for e in net.arcs
        for theta in 0:(T - e.transit_time)
            add_arc!(theta * n + e.tail, (theta + e.transit_time) * n + e.head, e.up)
        end
    end

    # 3. Super source/sink
    for theta in 0:T
        add_arc!(s_super, theta * n + s, c_rate)
        add_arc!(theta * n + t, t_super, c_rate)
    end

    # Build Sparse Matrix
    capacity_matrix = sparse(I, J, V, num_nodes_te, num_nodes_te)
    flow_graph = Graphs.SimpleDiGraph(capacity_matrix)

    # Run Dinic's Algorithm
    f_val, _ = GraphsFlows.maximum_flow(flow_graph, s_super, t_super, capacity_matrix, algorithm=GraphsFlows.DinicAlgorithm())
    
    return MDFResult(STATUS_OPTIMAL, Float64(f_val), Float64[], Tuple{Vector{Int},Float64}[])
end

"""
    solve_max_dynamic_flow_te_lp(net, s, t, T, solver, solver_method)

Solves the MDFP using a Time-Expanded network formulated as a Linear Program.

# Arguments
- `net`: The `NetworkInstance`.
- `s`: The source node.
- `t`: The sink node.
- `T`: The time horizon.
- `solver`: The JuMP solver backend (`"HiGHS"` or `"Gurobi"`).
- `solver_method`: The internal solver method (`"primal"`, `"dual"`, etc.).

# Returns
An `MDFResult`.
"""
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
    
    for (i, e) in enumerate(net.arcs)
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