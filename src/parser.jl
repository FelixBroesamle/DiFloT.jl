"""
    parse_min(filepath) -> NetworkInstance

Parse a DIMACS-style file. Robustly supports:
- `p min <nodes> <arcs> [target]`      (standard DIMACS from generators)
- `p <nodes> <arcs> [target]`           (some generators omit the keyword)
- `p MDFP <nodes> <arcs> <target>`      (DiFloT processed MDFP)
- `p QFP  <nodes> <arcs> <target>`      (DiFloT processed QFP)
"""
function parse_min(filepath::String)
    num_nodes = 0
    num_arcs = 0
    problem_type = :MIN
    target_value = 0.0
    source = 0
    sink = 0
    nodes = Int[]
    node_seen = Set{Int}()
    edges = Edge[]
    arcs = Tuple{Int,Int,Int}[]
    pair_counter = Dict{Tuple{Int,Int},Int}()
    arc_indices = Dict{Tuple{Int,Int},Vector{Int}}()

    open(filepath, "r") do f
        for raw in eachline(f)
            line = replace(String(raw), '\ufeff' => "")
            trimmed = strip(line)
            isempty(trimmed) && continue
            startswith(trimmed, 'c') && continue

            toks = split(trimmed)
            isempty(toks) && continue
            tag = uppercase(String(toks[1]))

            if tag == "P"
                length(toks) < 3 && error("Malformed p-line in $filepath: $line")

                # Detect whether toks[2] is a keyword or a number (old DIMACS)
                maybe_int = tryparse(Int, String(toks[2]))
                if maybe_int !== nothing
                    # Format: p <nodes> <arcs> [target]
                    problem_type = :MIN
                    num_nodes = maybe_int
                    num_arcs = parse(Int, String(toks[3]))
                    target_value = length(toks) >= 4 ? parse(Float64, String(toks[4])) : 0.0
                else
                    # Format: p <keyword> <nodes> <arcs> [target]
                    kw = lowercase(String(toks[2]))
                    length(toks) < 4 && error("Malformed p-line in $filepath: $line")
                    num_nodes = parse(Int, String(toks[3]))
                    num_arcs = parse(Int, String(toks[4]))
                    target_value = length(toks) >= 5 ? parse(Float64, String(toks[5])) : 0.0
                    if kw == "min"
                        problem_type = :MIN
                    elseif kw == "mdfp" || kw == "qfp"
                        problem_type = Symbol(uppercase(kw))
                    else
                        error("Unknown problem keyword '$kw' in $filepath")
                    end
                end
                sizehint!(edges, num_arcs)

            elseif tag == "N"
                nid = parse(Int, String(toks[2]))
                val = parse(Int, String(toks[3]))
                if val > 0 && source == 0; source = nid; end
                if val < 0 && sink == 0; sink = nid; end
                if !(nid in node_seen)
                    push!(node_seen, nid)
                    push!(nodes, nid)
                end

            elseif tag == "A"
                tail = parse(Int, String(toks[2]))
                head = parse(Int, String(toks[3]))
                low  = parse(Int, String(toks[4]))
                up   = parse(Int, String(toks[5]))
                tt   = parse(Int, String(toks[6]))

                slot = get(pair_counter, (tail, head), -1) + 1
                pair_counter[(tail, head)] = slot
                push!(edges, Edge(tail, head, low, up, tt, slot))

                for n in (tail, head)
                    if !(n in node_seen)
                        push!(node_seen, n)
                        push!(nodes, n)
                    end
                end

                push!(arcs, (tail, head, slot))
                push!(get!(arc_indices, (tail, head), Int[]), length(edges))
            end
        end
    end

    max_seen = isempty(nodes) ? 0 : maximum(nodes)
    for e in edges
        max_seen = max(max_seen, e.tail, e.head)
    end
    num_nodes = max(num_nodes, max_seen)

    return NetworkInstance(num_nodes, length(edges), sort!(nodes), edges,
                           problem_type, source, sink, target_value, arc_indices)
end

"""
    _format_target(v::Float64) -> String

Format the target value so integers are written without a decimal point
(e.g., `100` instead of `100.0`), keeping .mdfp / .qfp files idiomatic.
"""
function _format_target(v::Float64)
    if isinteger(v)
        return string(Int(v))
    else
        return string(v)
    end
end

"""
    save_min(net::NetworkInstance, filepath::String)

Save a `NetworkInstance` to disk in DIMACS-style format.

- If `net.problem_type == :MIN`, writes `p min <n> <m>` (no target).
- If `net.problem_type == :MDFP` / `:QFP`, writes
  `p MDFP <n> <m> <T>` / `p QFP <n> <m> <F>`.

Appends `.mdfp` or `.qfp` automatically if the filepath lacks the extension.
"""
function save_min(net::NetworkInstance, filepath::String)
    actual_path = filepath
    if net.problem_type != :MIN && !endswith(actual_path, ".mdfp") && !endswith(actual_path, ".qfp")
        actual_path *= (net.problem_type == :QFP ? ".qfp" : ".mdfp")
    end

    open(actual_path, "w") do io
        println(io, "c Generated by DiFloT.jl")
        if net.problem_type == :MIN
            println(io, "p min $(net.num_nodes) $(net.num_arcs)")
        else
            # FIX: always write target_value on processed MDFP/QFP instances
            println(io, "p $(net.problem_type) $(net.num_nodes) $(net.num_arcs) $(_format_target(net.target_value))")
        end

        if net.source > 0; println(io, "n $(net.source) 1");  end
        if net.sink   > 0; println(io, "n $(net.sink) -1");   end

        for e in net.edges
            println(io, "a $(e.tail) $(e.head) $(e.low) $(e.up) $(e.transit_time)")
        end
    end
end