"""
    solve_max_dynamic_flow_from_file(filepath; s=nothing, t=nothing, T=nothing, kwargs...)

Parse a DIMACS file and solve the Maximum Dynamic Flow Problem.
- `s`, `t`: If not provided, automatically extracts from the file.
- `T`: If not provided, reads the time horizon from a `.mdfp` file. Overrides the file value if provided.
"""
function solve_max_dynamic_flow_from_file(filepath::String; 
                                          s::Union{Nothing,Int}=nothing, 
                                          t::Union{Nothing,Int}=nothing, 
                                          T::Union{Nothing,Int}=nothing, 
                                          kwargs...)
    net = parse_min(filepath)
    
    actual_s = isnothing(s) ? net.source : s
    actual_t = isnothing(t) ? net.sink : t
    
    if actual_s == 0 || actual_t == 0
        error("Source/sink not found in file. Please provide `s` and `t` manually.")
    end

    actual_T = T
    if isnothing(actual_T)
        if net.problem_type == :MDFP
            actual_T = Int(net.target_value)
        else
            error("File is not an .mdfp instance (Found: $(net.problem_type)). Please provide a time horizon `T` explicitly.")
        end
    end
    
    return solve_max_dynamic_flow(net, actual_s, actual_t, actual_T; kwargs...)
end

"""
    solve_quickest_flow_from_file(filepath; s=nothing, t=nothing, F=nothing, kwargs...)

Parse a DIMACS file and solve the Quickest Flow Problem.
- `s`, `t`: If not provided, automatically extracts from the file.
- `F`: If not provided, reads the target flow from a `.qfp` file. Overrides the file value if provided.
"""
function solve_quickest_flow_from_file(filepath::String; 
                                       s::Union{Nothing,Int}=nothing, 
                                       t::Union{Nothing,Int}=nothing, 
                                       F::Union{Nothing,Number}=nothing, 
                                       kwargs...)
    net = parse_min(filepath)
    
    actual_s = isnothing(s) ? net.source : s
    actual_t = isnothing(t) ? net.sink : t
    
    if actual_s == 0 || actual_t == 0
        error("Source/sink not found in file. Please provide `s` and `t` manually.")
    end

    actual_F = F
    if isnothing(actual_F)
        if net.problem_type == :QFP
            actual_F = net.target_value
        else
            error("File is not a .qfp instance (Found: $(net.problem_type)). Please provide a flow target `F` explicitly.")
        end
    end
    
    return solve_quickest_flow(net, actual_s, actual_t, actual_F; kwargs...)
end