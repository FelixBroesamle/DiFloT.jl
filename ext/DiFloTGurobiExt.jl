module DiFloTGurobiExt

using Gurobi
import JuMP
import DiFloT: _gurobi_optimizer

function _gurobi_optimizer(sm::String)
    opt = Gurobi.Optimizer
    if sm == "primal"
        return JuMP.optimizer_with_attributes(opt, "Method" => 0)
    elseif sm == "dual"
        return JuMP.optimizer_with_attributes(opt, "Method" => 1)
    elseif sm == "network"
        return JuMP.optimizer_with_attributes(opt, "NetworkAlg" => 1)
    end
    return opt
end

end # module