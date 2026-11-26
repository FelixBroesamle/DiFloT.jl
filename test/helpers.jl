const METHODS = ("BS", "IBS", "INT", "INT-SG", "NEWT", "LP")

function repeated_value(net, T, paths)
    sum(delta * (T - sum(net.arcs[i].transit_time for i in P) + 1)
        for (P, delta) in paths; init=0.0)
end

function check_paths(net, s, t, T, paths)
    used = zeros(Float64, net.num_arcs)
    @test length(paths) <= net.num_arcs
    
    for (P, rate) in paths
        @test !isempty(P) && rate > 0
        node = s
        seen = Set([s])
        transit = 0
    
        for i in P
            e = net.arcs[i]
            @test e.tail == node
            @test !(e.head in seen)
            node = e.head
            push!(seen, node)
            transit += e.transit_time
            used[i] += rate
        end
    
        @test node == t && transit <= T
    end
    
    @test all(used[i] <= e.up + 1e-8 for (i,e) in enumerate(net.arcs))
end

function check_qfp(
    net, s, t, F, expected; 
    solver="HiGHS", solver_method="primal", max_flow_solver="GraphsFlows"
)
    for method in METHODS
        r = solve_quickest_flow(
            net, s, t, F; 
            method=method, solver=solver, solver_method=method == "LP" ? "dual" : solver_method, max_flow_solver=max_flow_solver,
            with_paths=true, verify=true
        )
        
        @test r.status == STATUS_OPTIMAL
        @test r.time_T == expected
        @test length(r.static_flow) == net.num_arcs
        @test repeated_value(net, Int(r.time_T), r.paths) >= F
        check_paths(net, s, t, Int(r.time_T), r.paths)
        @test r.mccp_calls >= 1
        @test r.lp_calls == (method == "LP" ? 1 : 0)
        @test r.safeguard_activations >= 0
    end
end

