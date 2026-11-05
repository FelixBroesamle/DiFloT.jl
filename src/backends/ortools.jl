module ORToolsBackend

using PythonCall
using ..DiFloT: NetworkInstance, Edge, MCCPResult, c_rate_bound, STATUS_OPTIMAL, STATUS_INFEASIBLE

const pygraph_mcf = PythonCall.pynew()
const pygraph_mf  = PythonCall.pynew()
const pynp        = PythonCall.pynew()
const _initialized = Ref(false)

function init_ortools!()
    if !_initialized[]
        PythonCall.pycopy!(pygraph_mcf, pyimport("ortools.graph.python.min_cost_flow"))
        PythonCall.pycopy!(pygraph_mf,  pyimport("ortools.graph.python.max_flow"))
        PythonCall.pycopy!(pynp,        pyimport("numpy"))
        _initialized[] = true
    end
end

function solve_max_flow(net, s::Int, t::Int)
    init_ortools!()
    starts = [e.tail - 1 for e in net.edges]
    ends   = [e.head - 1 for e in net.edges]
    caps   = [e.up       for e in net.edges]
    
    smf = pygraph_mf.SimpleMaxFlow()
    smf.add_arcs_with_capacity(pynp.array(starts), pynp.array(ends), pynp.array(caps))
    status = smf.solve(s - 1, t - 1)
    
    if pyconvert(Bool, status == smf.OPTIMAL)
        return (status = STATUS_OPTIMAL, value = pyconvert(Int64, smf.optimal_flow()))
    else
        return (status = STATUS_INFEASIBLE, value = 0)
    end
end

function solve_circulation(net, s::Int, t::Int, T::Int)
    init_ortools!()
    edges = copy(net.edges)
    c_rate = c_rate_bound(net, s, t)
    push!(edges, Edge(t, s, 0, c_rate, -(T + 1), 0)) 
    
    net_mod = NetworkInstance(net.num_nodes, length(edges), net.nodes, edges, net.problem_type, net.source, net.sink, net.target_value, net.arc_indices)
    
    smcf = pygraph_mcf.SimpleMinCostFlow()
    starts = [e.tail - 1 for e in net_mod.edges]
    ends   = [e.head - 1 for e in net_mod.edges]
    caps   = [e.up - e.low for e in net_mod.edges]
    costs  = [e.transit_time for e in net_mod.edges]
    
    supplies = zeros(Int, net_mod.num_nodes)
    base_cost = 0
    for (i, e) in enumerate(net_mod.edges)
        if e.low > 0
            supplies[e.tail] -= e.low
            supplies[e.head] += e.low
            base_cost += e.low * e.transit_time
        end
    end

    smcf.add_arcs_with_capacity_and_unit_cost(
        pynp.array(starts), pynp.array(ends), pynp.array(caps), pynp.array(costs))
    smcf.set_nodes_supplies(pynp.arange(0, net_mod.num_nodes), pynp.array(supplies))
    status = smcf.solve()

    if pyconvert(Bool, status == smcf.OPTIMAL)
        raw = pyconvert(Vector{Int64}, smcf.flows())
        flows = [raw[i] + net_mod.edges[i].low for i in 1:net_mod.num_arcs]
        v = -(pyconvert(Int64, smcf.optimal_cost()) + base_cost)
        g = flows[end]
        return MCCPResult(STATUS_OPTIMAL, v, g, flows[1:net.num_arcs])
    else
        return MCCPResult(STATUS_INFEASIBLE, -Inf, 0.0, Float64[])
    end
end

end # module