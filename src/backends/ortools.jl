# OR-Tools runs in a separate, persistent Python process. This prevents
# native-library conflicts between the OR-Tools wheel and HiGHS_jll.
const _ORTOOLS_LOCK = ReentrantLock()
const _ORTOOLS_IO = Ref{Union{Nothing, IO}}(nothing)
const _ORTOOLS_VERSION = Ref("")

"""
    start_ortools!(; python=get(ENV,"DIFLOT_PYTHON",Sys.iswindows() ? "python" : "python3"))

Start the optional OR-Tools worker and return its version. The Python
interpreter must have `ortools` and `numpy` installed. No packages are
installed automatically. Set DIFLOT_PYTHON to the interpreter path.
"""
function start_ortools!(; python::AbstractString=get(ENV, "DIFLOT_PYTHON", Sys.iswindows() ? "python" : "python3"))
    lock(_ORTOOLS_LOCK) do
        !isnothing(_ORTOOLS_IO[]) && return _ORTOOLS_VERSION[]
        Sys.WORD_SIZE == 64 && Base.ENDIAN_BOM == 0x04030201 ||
            error("OR-Tools IPC requires a little-endian 64-bit Julia process.")

        # Launch the worker and check its protocol handshake.
        script = joinpath(@__DIR__, "python", "ortools_worker.py")
        proc = open(pipeline(`$python -u $script`, stderr=stderr), "r+")
        try
            reply = readline(proc)
            startswith(reply, "DIFLOT-ORTOOLS/1 ") ||
                error("OR-Tools worker failed to initialize. Check DIFLOT_PYTHON and Python dependencies.")
            _ORTOOLS_VERSION[] = split(reply; limit=2)[2]
            _ORTOOLS_IO[] = proc
            return _ORTOOLS_VERSION[]

        catch
            close(proc)
            rethrow()
        end
    end
end

"""Stop the optional OR-Tools worker. A later request starts a fresh worker."""
function stop_ortools!()
    lock(_ORTOOLS_LOCK) do
        proc = _ORTOOLS_IO[]
        _ORTOOLS_IO[] = nothing
        _ORTOOLS_VERSION[] = ""
        if !isnothing(proc)
            try
                # Operation 0 requests a clean shutdown.
                write(proc, zeros(Int64, 5))
                flush(proc)
            catch
                # A failed worker may have already closed its pipes.
            finally
                close(proc)
            end
        end
    end
    return nothing
end

function _ortools_request(
    op::Int, I::Vector{Int}, J::Vector{Int}, V::Vector{Int}, n::Int, s::Int, t::Int; costs::Vector{Int}=Int[]
)
    m = length(V)
    length(I) == length(J) == m || error("Invalid arc arrays.")
    op == 1 || (op == 2 && length(costs) == m) || error("Invalid OR-Tools operation.")
    n <= typemax(Int32) || throw(ArgumentError("OR-Tools node index limit exceeded."))

    lock(_ORTOOLS_LOCK) do
        start_ortools!()
        proc = _ORTOOLS_IO[]::IO

        try
            # Send the header and arc arrays; Python uses zero-based node IDs.
            write(proc, Int64[op, n, m, s - 1, t - 1])
            write(proc, I .- 1)
            write(proc, J .- 1)
            write(proc, V)

            op == 2 && write(proc, costs)
            flush(proc)

            # Read the response header, followed by flows or an error message.
            status = read(proc, Int64)
            value = read(proc, Int64)
            count = read(proc, Int64)

            if status != 0
                0 <= count <= 1_000_000 || error("Invalid OR-Tools error response.")
                error("OR-Tools failed: " * String(read(proc, count)))
            end

            count == (op == 2 ? m : 0) || error("Invalid OR-Tools response length.")
            flows = read!(proc, Vector{Int64}(undef, count))

            return (value=Int(value), flows=flows)
        
        catch
            # Never reuse a potentially desynchronized stream.
            _ORTOOLS_IO[] = nothing
            _ORTOOLS_VERSION[] = ""

            try
                close(proc)
            catch
            end
            
            rethrow()
        end
    end
end

"""Solve the integer-horizon MCCP using OR-Tools."""
function _ortools_mccp(net, s::Int, t::Int, T::Int)
    _check_horizon(net, s, t, T)

    # Original arcs, in their original order, followed by the return arc.
    I = [e.tail for e in net.arcs]
    push!(I, t)
    J = [e.head for e in net.arcs]
    push!(J, s)
    V = [e.up for e in net.arcs]
    push!(V, c_rate_bound(net, s, t))
    costs = [e.transit_time for e in net.arcs]
    push!(costs, -(T + 1))

    res = _ortools_request(2, I, J, V, net.num_nodes, s, t; costs=costs)

    0 <= res.value <= MAX_EXACT_INT ||
        error("Invalid MCCP value.")

    return MCCPResult(
        STATUS_OPTIMAL,
        Float64(res.value),
        Float64(res.flows[end]),
        Float64.(res.flows[1:end-1])
    )
end
