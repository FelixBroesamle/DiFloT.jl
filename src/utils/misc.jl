using GraphsFlows

function _gurobi_optimizer(::String)
    error("Gurobi.jl is not loaded. Add and load it first:\n" *
          "    using Pkg; Pkg.add(\"Gurobi\")\n" *
          "    using Gurobi\n" *
          "then call solve_* with solver=\"Gurobi\".")
end


"""
    _optimizer(solver::String; solver_method::String="primal", warm_start::Bool=false)

Construct a JuMP optimizer with specific attributes based on the backend.
Supports `"HiGHS"` and `"Gurobi"`. Use `solver_method="primal"` for algorithms requiring warm-starts.
"""
function _optimizer(solver::String; solver_method::String="primal", warm_start::Bool=false)
    s = lowercase(solver)
    sm = lowercase(solver_method)

    if s == "highs"
        opt = HiGHS.Optimizer
        if warm_start || sm == "primal"
            # Primal simplex for warm-starting
            return JuMP.optimizer_with_attributes(opt, "solver" => "simplex", "simplex_strategy" => 1)
        elseif sm == "dual"
            return JuMP.optimizer_with_attributes(opt, "solver" => "simplex", "simplex_strategy" => 2)
        end
        return opt
    elseif s == "gurobi"
        return _gurobi_optimizer(sm) 
        #opt = Gurobi.Optimizer
        #if sm == "primal"
        #    return JuMP.optimizer_with_attributes(opt, "Method" => 0)
        #elseif sm == "dual"
        #    return JuMP.optimizer_with_attributes(opt, "Method" => 1)
        #elseif sm == "network"
        #    return JuMP.optimizer_with_attributes(opt, "NetworkAlg" => 1)
        #end
        #return opt
    end
    error("Unknown solver '$solver'.")
end

"""
    compute_max_static_flow(net, s, t)

Returns the maximum static flow value `g_max`.
"""
function compute_max_static_flow(net::NetworkInstance, s::Int, t::Int)
    n = net.num_nodes
    flow_graph = Graphs.SimpleDiGraph(n)
    capacity_matrix = zeros(Float64, n, n)
    
    for e in net.arcs
        if e.tail != e.head # Ignore self-loops
            Graphs.add_edge!(flow_graph, e.tail, e.head)
            capacity_matrix[e.tail, e.head] += e.up # Sum parallel capacities
        end
    end
    
    # GraphsFlows.maximum_flow returns (flow_value, flow_matrix)
    f_val, _ = GraphsFlows.maximum_flow(flow_graph, s, t, capacity_matrix)
    return (status = STATUS_OPTIMAL, value = f_val)
end
