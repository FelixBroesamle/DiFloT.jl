@testset "Example" begin
    net = NetworkInstance(
        [1,1,1,2,2,3,3,4], [2,3,4,3,5,5,4,5],
        [2,2,3,1,1,2,2,4], [1,3,7,1,5,3,1,1]; 
        source=1,sink=5
    )
    
    @test compute_shortest_path(net,1,5) == 4
    @test compute_max_static_flow(net,1,5).value == 7
    
    for (T, expected) in zip([0,3,4,5,6,7,8,9,15], [0,0,1,3,7,11,17,24,66])
        ff = solve_max_dynamic_flow(net,1,5,T;method="FF",solver_method="dual")
        te = solve_max_dynamic_flow(net,1,5,T;method="TE")
        @test ff.status == te.status == STATUS_OPTIMAL
        @test ff.value == te.value == expected
        @test repeated_value(net,T,ff.paths) == expected
        check_paths(net,1,5,T,ff.paths)
        @test all(sum(net.arcs[i].transit_time for i in P) <= T for (P,_) in ff.paths)
    end
    
    check_qfp(net,1,5,7,6)
    check_qfp(net,1,5,12,8)
    
    # Both simplex strategies should produce the same optimal value.
    @test solve_max_dynamic_flow(net,1,5,8;solver_method="primal").value == 17
    @test solve_maximum_dynamic_flow(net,1,5,8;with_paths=false).paths == []
    @test solve_quickest_flow(net,1,5,12;with_paths=false).paths == []
    
    lp = solve_quickest_flow(net,1,5,12;method="LP",verify=false)
    @test lp.lp_calls == 1 && lp.mccp_calls <= 2
end

@testset "Newton termination counterexample" begin
    # Four internally disjoint paths with capacities 88,32,128,217 and
    # transit times 1,2,4,7. The old step-length criterion returns 6.
    c = [88,32,128,217]
    net = NetworkInstance(
        [1,1,1,1,2,3,4,5],[2,3,4,5,6,6,6,6],
        vcat(c,c),[0,0,0,0,1,2,4,7];
        source=1,sink=6
    )
    @test solve_max_dynamic_flow(net,1,6,4).value == 576
    @test solve_max_dynamic_flow(net,1,6,5).value == 824
    check_qfp(net,1,6,700,5)
end

@testset "Zero transit times, parallel arcs and zero capacities" begin
    net = NetworkInstance([1,1,1,2],[2,2,3,3],[3,4,0,7],[0,0,0,0];source=1,sink=3)
    @test compute_shortest_path(net,1,3) == 0
    @test compute_max_static_flow(net,1,3).value == 7
    
    for T in 0:4
        @test solve_max_dynamic_flow(net,1,3,T;method="FF").value == 7*(T+1)
        @test solve_max_dynamic_flow(net,1,3,T;method="TE").value == 7*(T+1)
    end
    
    check_qfp(net,1,3,1,0)
    check_qfp(net,1,3,7,0)
    check_qfp(net,1,3,8,1)
    check_qfp(net,1,3,29,4)
    
    @test solve_quickest_flow(net,1,3,7;method="BS").mccp_calls == 1
    shortcut = NetworkInstance([1,1,2],[3,2,3],[0,2,2],[0,2,3];source=1,sink=3)
    @test compute_shortest_path(shortcut,1,3) == 5
    check_qfp(shortcut,1,3,3,6)
end

@testset "Independent TE cross-check on small networks" begin
    for k in 1:3
        net = NetworkInstance(
            [1,1,2,2,3,3,4,2],[2,3,3,4,2,4,5,5],
            [mod(i+k,5)+1 for i in 1:8],
            [mod(i*k,4) for i in 1:8];
            source=1,sink=5
        )
        
        values = [solve_max_dynamic_flow(net,1,5,T;method="TE").value for T in 0:10]
        
        for (T,v) in enumerate(values)
            @test solve_max_dynamic_flow(net,1,5,T-1;method="FF",with_paths=false).value == v
        end
        
        for F in (1,5,10)
            idx = findfirst(>=(F),values)
            if idx !== nothing
                check_qfp(net,1,5,F,idx-1)
            end
        end
    end
end

@testset "Safeguard activation" begin
    # Independent parallel paths; at least one tentative interpolation
    # step fails the 4/5 interval-reduction test for either extreme tie rule.
    c = [3,67,26,39,44,16,76]
    tt = [293,137,82,202,19,202,215]
    
    net = NetworkInstance(fill(1,7),fill(2,7),c,tt)
    
    r = solve_quickest_flow(net,1,2,12039;method="INT-SG",verify=true)
    
    @test r.time_T == 176
    @test r.safeguard_activations >= 1
end

@testset "LP-only and checked integer recovery" begin
    net = NetworkInstance([1],[2],[3],[1])
    F = 7
    c = solve_lp_candidate(net,1,2,F;solver_method="dual")
    
    @test c.status == STATUS_OPTIMAL && c.lp_calls == 1
    @test c.tau_tilde ≈ 10/3
    @test c.time_T == 3
    
    for T in (2,3,4) # Test both possible rounding errors, not just the usual case.
        tau = T == 2 ? 3.0 : T == 4 ? 4.1 : 10/3
        fake = LPCandidate(STATUS_OPTIMAL,tau,T,1)
        r = check_lp_candidate(net,1,2,F,fake;solver_method="primal")
        @test r.time_T == 3 && r.mccp_calls <= 2 && r.lp_calls == 1
    end
    
    @test solve_quickest_flow(net,1,2,F;method="LP-check",solver_method="dual").time_T == 3
end
