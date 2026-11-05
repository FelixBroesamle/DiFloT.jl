"""
    solve_quickest_flow_from_file(filepath, s, t, F; kwargs...)

Parse a DIMACS `.min` file and solve the Quickest Flow Problem directly.
"""
function solve_quickest_flow_from_file(filepath::String, s::Int, t::Int, F::Number; kwargs...)
    net = parse_min(filepath)
    return solve_quickest_flow(net, s, t, F; kwargs...)
end