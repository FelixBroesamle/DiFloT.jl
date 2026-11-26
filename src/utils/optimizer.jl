"""Load an optional backend without overwriting a method during precompilation."""
function _require_extension(name::Symbol, package::String)
    ext = Base.get_extension(@__MODULE__, name)
    isnothing(ext) && error("Load $package first (using $package) in an environment where it is installed.")
    return ext
end

_gurobi_optimizer(sm::String) = _require_extension(:DiFloTGurobiExt, "Gurobi").make_optimizer(sm)

"""
    _optimizer(solver; solver_method="primal", warm_start=false)

Construct a single-threaded optimizer. Reuse of an MCCP model is independent
of the chosen simplex strategy; basis reuse is backend-dependent.
"""
function _optimizer(solver::String; solver_method::String="primal", warm_start::Bool=false)
    s, sm = lowercase(solver), lowercase(solver_method)
    sm in ("primal", "dual", "network", "auto") ||
        throw(ArgumentError("Unknown solver method '$solver_method'."))
    if s == "highs"
        sm == "network" && throw(ArgumentError("HiGHS has no network simplex option."))
        strategy = sm == "primal" ? 4 : (sm == "dual" ? 1 : 0)
        return JuMP.optimizer_with_attributes(HiGHS.Optimizer,
                    "solver" => "simplex", "simplex_strategy" => strategy,
                    "threads" => 1)
    elseif s == "gurobi"
        return _gurobi_optimizer(sm)
    end
    throw(ArgumentError("Unknown solver '$solver'."))
end
