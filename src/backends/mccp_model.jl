# Integer oracle backed by the persistent OR-Tools worker.
mutable struct ORToolsMCCPSession
    net::NetworkInstance
    s::Int
    t::Int
    last_T::Int
    last_result::Union{Nothing,MCCPResult}
    calls::Int
end

function solve_mccp(session::ORToolsMCCPSession, theta::Real; force::Bool=false)
    T = _exact_int(theta, "OR-Tools MCCP parameter")
    _check_horizon(session.net, session.s, session.t, T)

    # Reuse the last evaluation without counting a new oracle call. 
    if !force && session.last_T == T && session.last_result !== nothing
        return session.last_result
    end

    # Invalidate the cache before solving.
    session.last_result = nothing
    session.calls += 1
    r = _ortools_mccp(session.net, session.s, session.t, T)
    session.last_T = T
    session.last_result = r
    return r
end

# Reusable JuMP model for integer and real-valued MCCP parameters.
mutable struct MCCPSession
    net::NetworkInstance
    s::Int
    t::Int
    model::JuMP.Model
    x::Vector{JuMP.VariableRef}
    x_ret::JuMP.VariableRef
    last_coefficient::Float64
    last_theta::Float64
    last_result::Union{Nothing,MCCPResult}
    calls::Int
    rate_bound::Int
end

"""
    MCCPSession(net, s, t; solver="HiGHS", solver_method="primal")

Create an MCCP oracle session. LP backends reuse the model; OR-Tools
solves a new integer circulation at each uncached evaluation.
"""
function MCCPSession(
    net::NetworkInstance,s::Int,t::Int;
    solver::String="HiGHS",solver_method::String="primal"
)
    _validate_network(net,s,t)
    if lowercase(solver) == "ortools"
        return ORToolsMCCPSession(net, s, t, -1, nothing, 0)
    end

    # Build the model once; subsequent calls change only the return-arc cost.
    model = Model(_optimizer(solver;solver_method=solver_method,warm_start=true))
    set_silent(model)

    # One variable per original arc, preserving parallel arcs.
    x = JuMP.VariableRef[]
    for e in net.arcs
        push!(x,@variable(model,lower_bound=0.0,upper_bound=Float64(e.up)))
    end

    # Add the return arc (t, s).
    c_rate = c_rate_bound(net,s,t)
    x_ret = @variable(model,lower_bound=0.0,upper_bound=Float64(c_rate))

    # Flow conservation on the augmented network.
    incoming,outgoing = positional_incidence(net)
    for v in 1:net.num_nodes
        inflow = sum(x[i] for i in incoming[v];init=0.0)
        outflow = sum(x[i] for i in outgoing[v];init=0.0)
        if v == s; inflow += x_ret; end
        if v == t; outflow += x_ret; end
        @constraint(model,inflow == outflow)
    end

    @objective(model,Min,sum(net.arcs[i].transit_time*x[i] for i in eachindex(x);init=0.0))
    
    return MCCPSession(net,s,t,model,x,x_ret,NaN,NaN,nothing,0,c_rate)
end

"""
    solve_mccp(session, theta; tol=FLOW_TOL, force=false) -> MCCPResult

Evaluate the continuous value function at `theta`. For integer parameters,
reconstruct the integer objective from a rounded, exactly feasible static
flow. A nonoptimal solver status or invalid flow is an error, not evidence
of network infeasibility. Calls count actual optimize! invocations.

A numerical LP optimum is not an independent exact optimality certificate.
"""
function solve_mccp(session::MCCPSession,theta::Real;tol::Real=FLOW_TOL,force::Bool=false)
    isfinite(theta) && 0 <= theta < MAX_EXACT_INT || throw(ArgumentError("Invalid MCCP parameter $theta."))
    if isinteger(theta)
        Int128(Int(theta)+1)*session.rate_bound <= MAX_EXACT_INT ||
            throw(ArgumentError("Horizon/capacity product exceeds exact-integer range."))
    else
        (Float64(theta)+1)*session.rate_bound <= MAX_EXACT_INT ||
            throw(ArgumentError("MCCP parameter exceeds supported range."))
    end

    # The only parameter-dependent cost is -(theta + 1) on the return arc.
    a = -(Float64(theta)+1.0)
    if !force && session.last_theta == Float64(theta) && session.last_result !== nothing
        return session.last_result
    end

    # Invalidate the old result before attempting a new solve.
    session.last_coefficient = NaN
    session.last_theta = NaN
    session.last_result = nothing

    # Reoptimize the existing model.
    JuMP.set_objective_coefficient(session.model,session.x_ret,a)
    session.calls += 1
    optimize!(session.model)
    _require_optimal(session.model,"MCCP($theta)")

    # Recover integral arc flows and verify capacities and balances.
    flows = Vector{Float64}(undef,session.net.num_arcs)
    balance = zeros(Int,session.net.num_nodes)
    cost = Int128(0)
    for (i,e) in enumerate(session.net.arcs)
        xi = _rounded_flow(value(session.x[i]),"MCCP arc flow";tol=tol)
        0 <= xi <= e.up || error("MCCP flow violates capacity on arc $i.")
        flows[i] = Float64(xi)

        balance[e.tail] += xi
        balance[e.head] -= xi

        cost += Int128(e.transit_time)*xi
    end

    # The return-arc flow is both the static flow value and a subgradient.
    g = _rounded_flow(value(session.x_ret),"return-arc flow";tol=tol)
    0 <= g <= session.rate_bound || error("Invalid MCCP return-arc flow.")

    for v in eachindex(balance)
        balance[v] == (v == session.s ? g : v == session.t ? -g : 0) ||
            error("Rounded MCCP flow violates conservation at node $v.")
    end

    # Reconstruct v(theta) = (theta + 1)*g - cost.
    v = if isinteger(theta)
        z = (Int128(Int(theta))+1)*g-cost
        0 <= z <= MAX_EXACT_INT || error("Invalid integer MCCP value $z.")
        Float64(z)
    else
        z = muladd(Float64(theta),Float64(g),Float64(Int128(g)-cost))
        z >= -tol || error("Negative MCCP value $z.")
        max(0.0,z)
    end

    # Cache only the successfully evaluated result.
    session.last_coefficient = a
    session.last_theta = Float64(theta)
    session.last_result = MCCPResult(STATUS_OPTIMAL,v,Float64(g),flows)
    
    return session.last_result
end
