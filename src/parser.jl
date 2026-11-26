"""
    parse_min(filepath) -> NetworkInstance

Parse a DIMACS-style file. Supports:
- `p min <nodes> <arcs> [target]`
- `p <nodes> <arcs> [target]`
- `p MDFP <nodes> <arcs> <target>`
- `p QFP <nodes> <arcs> <target>`

Preserves arc ordering, orientation, parallel arcs and isolated nodes.
"""
function parse_min(filepath::String)
    num_nodes = 0
    num_arcs = 0
    problem_type = :MIN
    target_value = 0.0
    source = 0
    sink = 0
    arcs = Arc[]
    arc_indices = Dict{Tuple{Int,Int},Vector{Int}}()
    seen_header = false
    node_values = Dict{Int,Int}()

    open(filepath, "r") do f
        for (lineno, raw) in enumerate(eachline(f))
            line = strip(replace(raw, '\ufeff' => ""))
            isempty(line) && continue
            toks = split(line)
            tag = uppercase(toks[1])
            tag == "C" && continue
            if tag == "P"
                seen_header && error("Duplicate p-line in $filepath:$lineno.")
                length(toks) >= 3 || error("Malformed p-line in $filepath:$lineno.")
                maybe_n = tryparse(Int, toks[2])
                if isnothing(maybe_n)
                    length(toks) >= 4 || error("Malformed p-line in $filepath:$lineno.")
                    problem_type = Symbol(uppercase(toks[2]))
                    problem_type in (:MIN, :MDFP, :QFP) || error("Unknown problem type $problem_type.")
                    num_nodes = parse(Int, toks[3])
                    num_arcs = parse(Int, toks[4])
                    target_idx = 5
                else
                    num_nodes = maybe_n
                    num_arcs = parse(Int, toks[3])
                    target_idx = 4
                end
                num_nodes > 0 && num_arcs >= 0 || error("Invalid network size.")
                length(toks) in (target_idx-1, target_idx) ||
                    error("Malformed p-line in $filepath:$lineno.")
                if length(toks) >= target_idx
                    target = parse(BigFloat, toks[target_idx])
                    target_value = problem_type == :MIN ? _target_float(target) :
                                   Float64(_exact_int(target, "target"))
                elseif problem_type != :MIN
                    error("Missing target in $filepath:$lineno.")
                end
                if problem_type == :MDFP
                    target_value >= 0 || error("MDFP requires a nonnegative horizon.")
                elseif problem_type == :QFP
                    target_value > 0 || error("QFP requires a positive demand.")
                end
                sizehint!(arcs, num_arcs)
                seen_header = true
            elseif tag == "N"
                seen_header && length(toks) == 3 || error("Malformed n-line in $filepath:$lineno.")
                nid, val = parse(Int, toks[2]), parse(Int, toks[3])
                1 <= nid <= num_nodes || error("Node ID out of range in $filepath:$lineno.")
                get(node_values, nid, val) == val ||
                    error("Conflicting node declarations in $filepath:$lineno.")
                node_values[nid] = val
                if val > 0
                    source == 0 || source == nid || error("Multiple sources are not supported.")
                    source = nid
                elseif val < 0
                    sink == 0 || sink == nid || error("Multiple sinks are not supported.")
                    sink = nid
                end
            elseif tag == "A"
                seen_header && length(toks) == 6 || error("Malformed a-line in $filepath:$lineno.")
                tail, head, low, up, tt = parse.(Int, toks[2:6])
                1 <= tail <= num_nodes && 1 <= head <= num_nodes || error("Invalid arc endpoint in $filepath:$lineno.")
                0 <= low <= up && tt >= 0 || error("Invalid arc data in $filepath:$lineno.")
                indices = get!(arc_indices, (tail, head), Int[])
                push!(arcs, Arc(tail, head, low, up, tt, length(indices)))
                push!(indices, length(arcs))
            else
                error("Unknown record '$tag' in $filepath:$lineno.")
            end
        end
    end
    seen_header || error("Missing p-line in $filepath.")
    length(arcs) == num_arcs || error("Declared and actual arc counts differ in $filepath.")
    return NetworkInstance(num_nodes, num_arcs, collect(1:num_nodes), arcs,
                           problem_type, source, sink, target_value, arc_indices)
end

_format_target(v::Float64) = isinteger(v) && abs(v) <= MAX_EXACT_INT ? string(Int(v)) : string(v)

"""
    save_min(net, filepath)

Write a DIMACS-style network file.
"""
function save_min(net::NetworkInstance, filepath::String)
    net.problem_type in (:MIN, :MDFP, :QFP) || error("Invalid problem type.")
    expected_ext = net.problem_type == :MIN ? ".min" :
                   net.problem_type == :MDFP ? ".mdfp" : ".qfp"
    ext = lowercase(splitext(filepath)[2])
    isempty(ext) || ext == expected_ext || error("Expected extension '$expected_ext'.")
    actual_path = isempty(ext) ? filepath * expected_ext : filepath
    net.num_arcs == length(net.arcs) || error("Inconsistent arc count.")
    net.num_nodes > 0 || error("Invalid node count.")
    if net.problem_type == :MDFP
        _exact_int(net.target_value, "T") >= 0 || error("Invalid time-horizon.")
    elseif net.problem_type == :QFP
        _positive_demand(net.target_value)
    end
    for e in net.arcs
        1 <= e.tail <= net.num_nodes && 1 <= e.head <= net.num_nodes || error("Invalid arc endpoint.")
        0 <= e.low <= e.up && e.transit_time >= 0 || error("Invalid arc data.")
    end
    0 <= net.source <= net.num_nodes && 0 <= net.sink <= net.num_nodes || error("Invalid terminal ID.")
    (net.source == 0 || net.sink == 0 || net.source != net.sink) || error("Source and sink must be distinct.")
    open(actual_path, "w") do io
        println(io, "c Generated by DiFloT.jl")
        if net.problem_type == :MIN
            println(io, "p min $(net.num_nodes) $(net.num_arcs)")
        else
            println(io, "p $(net.problem_type) $(net.num_nodes) $(net.num_arcs) $(_format_target(net.target_value))")
        end
        net.source > 0 && println(io, "n $(net.source) 1")
        net.sink > 0 && println(io, "n $(net.sink) -1")
        for e in net.arcs
            println(io, "a $(e.tail) $(e.head) $(e.low) $(e.up) $(e.transit_time)")
        end
    end
    return nothing
end
