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

Returns the maximum static flow value `g_max` using HiGHS (no Python required).
"""
function compute_max_static_flow(net, s::Int, t::Int)
    model = Model(_optimizer("HiGHS", solver_method="primal"))
    set_silent(model)
    
    x = Dict{Tuple{Int,Int,Int},JuMP.VariableRef}()
    for e in net.edges
        k = arc_key(e)
        x[k] = @variable(model, lower_bound=Float64(e.low), upper_bound=Float64(e.up))
    end
    
    in_keys, out_keys = incidence_mapping(net)
    for v in 1:net.num_nodes
        if v == s || v == t; continue; end
        @constraint(model, sum(x[k] for k in in_keys[v]; init=0.0) == sum(x[k] for k in out_keys[v]; init=0.0))
    end
    
    @objective(model, Max, sum(x[k] for k in out_keys[s]; init=0.0))
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        return (status = STATUS_OPTIMAL, value = objective_value(model))
    else
        return (status = STATUS_INFEASIBLE, value = 0.0)
    end
end