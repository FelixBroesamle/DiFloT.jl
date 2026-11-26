"""
    solve_lp_candidate(net,s,t,F; solver="HiGHS", solver_method="dual")

Solve the continuous Lin--Jaillet LP and round its objective to a discrete
horizon. No MCCP is solved. STATUS_OPTIMAL refers only to the numerical LP
solve: the discrete candidate is not verified. No static max-flow precheck
is needed because the LP detects infeasibility directly.
"""
function solve_lp_candidate(net,s::Int,t::Int,F::Real;solver::String="HiGHS",
                            solver_method::String="dual")
    F_int = _positive_demand(F)
    _validate_network(net,s,t)
    lowercase(solver_method) == "network" && throw(ArgumentError("The Lin-Jaillet LP is not a network simplex model."))
    model = Model(_optimizer(solver;solver_method=solver_method))
    set_silent(model)
    x = JuMP.VariableRef[]
    for e in net.arcs
        push!(x,@variable(model,lower_bound=0.0))
    end
    theta = @variable(model,lower_bound=0.0)
    @objective(model,Min,F_int*theta+sum(net.arcs[i].transit_time*x[i] for i in eachindex(x);init=0.0))
    incoming,outgoing = positional_incidence(net)
    for v in 1:net.num_nodes
        rhs = v == s ? 1.0 : v == t ? -1.0 : 0.0
        @constraint(model,sum(x[i] for i in outgoing[v];init=0.0)-
                          sum(x[i] for i in incoming[v];init=0.0) == rhs)
    end
    for (i,e) in enumerate(net.arcs)
        @constraint(model,x[i] <= e.up*theta)
    end
    optimize!(model)
    if termination_status(model) == MOI.INFEASIBLE
        return LPCandidate(STATUS_INFEASIBLE,Inf,-1,1)
    end
    _require_optimal(model,"Lin-Jaillet LP")
    tau_tilde = objective_value(model)
    isfinite(tau_tilde) && -1 < tau_tilde <= MAX_EXACT_INT || error("Invalid LP objective $tau_tilde.")
    T = max(0,ceil(Int,tau_tilde)-1)

    return LPCandidate(STATUS_OPTIMAL,tau_tilde,T,1)
end

"""
    check_lp_candidate(net,s,t,F,candidate; solver="HiGHS", solver_method="primal", with_paths=false, verify=false)

Recover the discrete optimum using at most two MCCP solves under the paper's
LP accuracy assumption. The candidate must come from the same instance.
With verify=true also check the predecessor if it was not already checked.
"""
function check_lp_candidate(net,s::Int,t::Int,F::Real,candidate::LPCandidate;
                            solver::String="HiGHS",solver_method::String="primal",
                            with_paths::Bool=false,verify::Bool=false)
    F_int = _positive_demand(F)
    _validate_network(net,s,t)
    if candidate.status == STATUS_INFEASIBLE
        return QFResult(STATUS_INFEASIBLE,-1.0,Float64[],PathFlow[],0,
                        candidate.lp_calls,0)
    end
    candidate.status == STATUS_OPTIMAL || error("Invalid LP candidate status.")
    T = candidate.time_T
    _check_horizon(net,s,t,T)
    session = MCCPSession(net,s,t;solver=solver,solver_method=solver_method)
    r = solve_mccp(session,T)
    predecessor_checked = false
    if _value(r) < F_int
        predecessor_checked = true
        T += 1
        r = solve_mccp(session,T)
    elseif T > 0
        r_prev = solve_mccp(session,T-1)
        if _value(r_prev) >= F_int
            T -= 1
            r = r_prev
        else
            predecessor_checked = true
        end
    end
    return _qf_result(session,T,r,F_int;with_paths=with_paths,verify=verify,lp_calls=candidate.lp_calls,predecessor_checked=predecessor_checked)
end

"""Lin--Jaillet LP followed by MCCP recovery (LP-check)."""
function solve_qf_lp(net,s::Int,t::Int,F::Real;solver::String="HiGHS",
                     solver_method::String="dual",max_flow_solver::String="GraphsFlows",
                     mccp_solver::String=solver,mccp_solver_method::String="primal",
                     with_paths::Bool=false,verify::Bool=false)
    candidate = solve_lp_candidate(net,s,t,F;solver=solver,solver_method=solver_method)
    return check_lp_candidate(net,s,t,F,candidate;solver=mccp_solver,
                              solver_method=mccp_solver_method,with_paths=with_paths,verify=verify)
end
