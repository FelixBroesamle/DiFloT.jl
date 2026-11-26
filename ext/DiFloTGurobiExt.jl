module DiFloTGurobiExt

using Gurobi
import JuMP

function make_optimizer(sm::String)
    opt = Gurobi.Optimizer
    if sm == "primal"
        return JuMP.optimizer_with_attributes(opt, "Method" => 0, "Threads" => 1)
    elseif sm == "dual"
        return JuMP.optimizer_with_attributes(opt, "Method" => 1, "Threads" => 1)
    elseif sm == "network"
        throw(ArgumentError("Network simplex is not supported by this backend configuration; choose primal or dual."))
    end
    return JuMP.optimizer_with_attributes(opt, "Method" => 1, "Threads" => 1)
end

end # module
