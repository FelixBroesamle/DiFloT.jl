"""
    solve_max_dynamic_flow_from_file(filepath; s=nothing, t=nothing, T=nothing, kwargs...)

Parse a DIMACS-style file and solve the MDFP. Explicit arguments override
source, sink, and time-horizon values read from the file.
"""
function solve_max_dynamic_flow_from_file(filepath::String;
                                          s::Union{Nothing,Int}=nothing,
                                          t::Union{Nothing,Int}=nothing,
                                          T::Union{Nothing,Int}=nothing,
                                          kwargs...)
    net = parse_min(filepath)
    actual_s = isnothing(s) ? net.source : s
    actual_t = isnothing(t) ? net.sink : t
    actual_s != 0 && actual_t != 0 || error("Provide source and sink explicitly.")
    if isnothing(T)
        net.problem_type == :MDFP || error("Provide a time-horizon T explicitly.")
        T = _exact_int(net.target_value, "T")
    end
    return solve_max_dynamic_flow(net, actual_s, actual_t, T; kwargs...)
end

"""
    solve_quickest_flow_from_file(filepath; s=nothing, t=nothing, F=nothing, kwargs...)

Parse a DIMACS-style file and solve the QFP. Explicit arguments override
source, sink, and demand values read from the file.
"""
function solve_quickest_flow_from_file(filepath::String;
                                       s::Union{Nothing,Int}=nothing,
                                       t::Union{Nothing,Int}=nothing,
                                       F::Union{Nothing,Real}=nothing,
                                       kwargs...)
    net = parse_min(filepath)
    actual_s = isnothing(s) ? net.source : s
    actual_t = isnothing(t) ? net.sink : t
    actual_s != 0 && actual_t != 0 || error("Provide source and sink explicitly.")
    if isnothing(F)
        net.problem_type == :QFP || error("Provide a demand F explicitly.")
        F = net.target_value
    end
    return solve_quickest_flow(net, actual_s, actual_t, F; kwargs...)
end
