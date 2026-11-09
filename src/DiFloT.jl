module DiFloT

using Printf
using Graphs
import HiGHS
#import Gurobi
import JuMP
using JuMP: Model, VariableRef, @variable, @constraint, @objective,
            optimize!, termination_status, objective_value, value,
            set_silent, set_optimizer_attribute, MOI

# Core definitions and types
include("types.jl")
include("parser.jl")

# Utilities
include("utils/misc.jl")
include("utils/mappings.jl")
include("utils/capacity_bounds.jl")
include("utils/shortest_path.jl")
include("utils/path_decomposition.jl")

# Backends
include("backends/mccp_model.jl")

# MDF Algorithms
include("mdfp/mdf.jl")
include("mdfp/time_expanded.jl")

# QF Algorithms
include("qfp/qf.jl")
include("qfp/binary_search.jl")
include("qfp/improved_binary_search.jl")
include("qfp/interpolation_search.jl")
include("qfp/newton.jl")
include("qfp/lp.jl")

# Direct File interface
include("direct.jl")

# Exports
export Arc, NetworkInstance
export parse_min, save_min
export solve_max_dynamic_flow, solve_quickest_flow
export solve_quickest_flow_from_file, solve_max_dynamic_flow_from_file
export MDFResult, QFResult, STATUS_OPTIMAL, STATUS_INFEASIBLE
export compute_shortest_path, compute_max_static_flow

end # module
