"""
    _gurobi_optimizer(::String)

Internal fallback that throws an error if Gurobi is requested but the extension is not loaded.
"""
function _gurobi_optimizer(::String)
    error("Gurobi.jl is not loaded. Add and load it first:\n" *
          "    using Pkg; Pkg.add(\"Gurobi\")\n" *
          "    using Gurobi\n" *
          "then call solve_* with solver=\"Gurobi\".")
end

"""
    _optimizer(solver::String; solver_method::String="primal", warm_start::Bool=false)

Construct a JuMP optimizer with specific attributes based on the backend.
Supports `"HiGHS"` (default) and `"Gurobi"` (via the package extension).
"""
function _optimizer(solver::String; solver_method::String="primal", warm_start::Bool=false)
    s = lowercase(solver)
    sm = lowercase(solver_method)

    if s == "highs"
        opt = HiGHS.Optimizer
        if warm_start || sm == "primal"
            return JuMP.optimizer_with_attributes(opt, "solver" => "simplex", "simplex_strategy" => 1)
        elseif sm == "dual"
            return JuMP.optimizer_with_attributes(opt, "solver" => "simplex", "simplex_strategy" => 2)
        end
        return opt
    elseif s == "gurobi"
        # This dispatches to the extension if Gurobi is loaded.
        # If Gurobi is not loaded, it triggers the fallback error above.
        return _gurobi_optimizer(sm) 
    end
    error("Unknown solver '$solver'.")
end