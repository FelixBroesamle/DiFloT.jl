module DiFloT

using Printf
using Graphs
import GraphsFlows
using SparseArrays
import HiGHS
import JuMP
using JuMP: Model, VariableRef, @variable, @constraint, @objective,
            optimize!, termination_status, objective_value, value,
            set_silent, set_optimizer_attribute, MOI

# Core definitions and types
include("types.jl")
include("parser.jl")

# Utilities
include("utils/validation.jl")
include("utils/optimizer.jl")
include("backends/ortools.jl")
include("utils/max_flow.jl")
include("utils/mappings.jl")
include("utils/capacity_bounds.jl")
include("utils/shortest_path.jl")
include("utils/path_decomposition.jl")

# Backends
include("backends/mccp_model.jl")

# MDFP algorithms
include("mdfp/mdf.jl")
include("mdfp/time_expanded.jl")

# QFP algorithms
include("qfp/common.jl")
include("qfp/qf.jl")
include("qfp/binary_search.jl")
include("qfp/improved_binary_search.jl")
include("qfp/interpolation_search.jl")
include("qfp/newton.jl")
include("qfp/lp.jl")

# Direct file interface
include("direct.jl")

export Arc, NetworkInstance
export parse_min, save_min
export solve_max_dynamic_flow, solve_maximum_dynamic_flow, solve_quickest_flow
export solve_quickest_flow_from_file, solve_max_dynamic_flow_from_file
export MDFResult, QFResult, STATUS_OPTIMAL, STATUS_INFEASIBLE
export compute_shortest_path, compute_max_static_flow
export solve_lp_candidate, check_lp_candidate, LPCandidate
export start_ortools!, stop_ortools!

function __init__()
    atexit(stop_ortools!)
end

end # module
