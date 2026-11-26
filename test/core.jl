@testset "Types and constructor" begin
    @test STATUS_OPTIMAL == 0
    @test STATUS_INFEASIBLE == 1
    
    @test fieldnames(MDFResult) == (:status, :value, :flows, :paths)
    @test fieldnames(QFResult) == (:status, :time_T, :static_flow, :paths, :mccp_calls, :lp_calls, :safeguard_activations)
    @test fieldnames(LPCandidate) == (:status, :tau_tilde, :time_T, :lp_calls)
    
    e = Arc(1, 2, 0, 10, 1)
    
    @test (e.tail, e.head, e.low, e.up, e.transit_time, e.parallel_index) == (1,2,0,10,1,0)
    
    net = NetworkInstance([1,1,2,3],[2,3,4,4],[10,5,10,10],[1,2,1,1];source=1,sink=4)
    
    @test net.num_nodes == 4 && net.num_arcs == length(net.arcs) == 4
    @test net.arc_indices[(1,2)] == [1]
    @test compute_shortest_path(net,1,4) == 2
    @test compute_shortest_path(net,1,1) == 0
end

@testset "Infeasible networks and invalid inputs" begin
    net = NetworkInstance([1],[2],[1],[1];source=1,sink=3)
    
    @test compute_shortest_path(net,1,3) == Inf
    @test compute_max_static_flow(net,1,3).value == 0
    
    for method in METHODS
        @test solve_quickest_flow(net,1,3,1;method=method).status == STATUS_INFEASIBLE
        for F in (0,-1,1.5,Inf,NaN)
            @test_throws ArgumentError solve_quickest_flow(net,1,3,F;method=method)
        end
    end
    
    @test solve_max_dynamic_flow(net,1,3,20).value == 0
    @test_throws ArgumentError solve_max_dynamic_flow(net,1,3,-1)
    @test_throws ArgumentError solve_quickest_flow(net,1,1,1)
    @test_throws ArgumentError solve_quickest_flow(net,1,4,1)
    @test_throws ArgumentError NetworkInstance([1],[2],[-1],[1])
end

@testset "Path decomposition and filtering" begin
    # Positive net flow 1 plus a circulation of value 1. A naive extraction
    # using only the path bottleneck would send 2 units instead of 1.
    net = NetworkInstance([1,2],[2,1],[2,1],[0,0];source=1,sink=2)
    paths = DiFloT.decompose_flow_into_paths(net,1,2,[2.0,1.0])
    @test paths == [([1],1.0)]
    paths = DiFloT.decompose_flow_into_paths(net,1,2,[1.5,1.0])
    @test paths == [([1],0.5)]
    @test isempty(DiFloT.decompose_flow_into_paths(net,1,2,[1.0,1.0]))
    
    # At T=0, transit time 1 has zero contribution to the MCCP optimum.
    net2 = NetworkInstance([1,1],[2,2],[1,1],[0,1];source=1,sink=2)
    @test DiFloT.decompose_flow_into_paths(net2,1,2,[1.0,1.0];T=0) == [([1],1.0)]
end

@testset "DIMACS round trip" begin
    mktempdir() do dir
        p = joinpath(dir,"example.qfp")
        write(p,"c test\np QFP 4 2 7\nn 1 1\nn 2 -1\na 1 2 0 3 0\na 1 2 0 4 1\n")
    
        net = parse_min(p)
        @test net.nodes == [1,2,3,4]
        @test net.arcs[2].parallel_index == 1
        @test solve_quickest_flow_from_file(p).time_T == 1
    
        q = joinpath(dir,"copy.qfp")
        save_min(net,q)
    
        other = parse_min(q)
    
        @test other.target_value == net.target_value
        @test other.arc_indices == net.arc_indices
        write(p,"p QFP 2 1 1\nn 1 1\nn 2 -1\na 1 2 1 3 0\n")
        @test_throws ArgumentError solve_quickest_flow_from_file(p)
        write(p,"p QFP 2 2 1\nn 1 1\nn 2 -1\na 1 2 0 3 0\n")
        @test_throws ErrorException parse_min(p)
    end
end

@testset "Directed and parallel arcs" begin
    mktempdir() do dir
        file = joinpath(dir, "parallel.min")
        write(file, """
        p min 4 4
        n 1 1
        n 3 -1
        a 1 2 0 3 0
        a 2 1 0 99 7
        a 1 2 0 5 2
        a 2 3 0 8 1
        """)
        net = parse_min(file)
        
        arc_data(n) = [(a.tail, a.head, a.low, a.up, a.transit_time, a.parallel_index) for a in n.arcs]
        
        @test net.nodes == [1, 2, 3, 4]
        @test arc_data(net) == [(1,2,0,3,0,0), (2,1,0,99,7,0), (1,2,0,5,2,1), (2,3,0,8,1,0)]
        @test net.arc_indices == Dict((1,2)=>[1,3], (2,1)=>[2], (2,3)=>[4])
        
        incoming, outgoing = DiFloT.positional_incidence(net)
        @test incoming == [[2], [1,3], [4], Int[]]
        @test outgoing == [[1,3], [2,4], Int[], Int[]]
        
        out = joinpath(dir, "copy.min")
        
        save_min(net, out)
        
        @test arc_data(parse_min(out)) == arc_data(net)
        @test compute_shortest_path(net, 1, 3) == 1
        @test compute_shortest_path(net, 3, 1) == Inf
        @test compute_max_static_flow(net, 1, 3).value == 8
        
        for T in 0:5
            expected = 3*T + 5*max(0,T-2)
            @test solve_max_dynamic_flow(net, 1, 3, T; method="FF").value == expected
            @test solve_max_dynamic_flow(net, 1, 3, T; method="TE").value == expected
        end
        
        write(file, "p QFP 2 0 0\n")
        
        @test_throws ErrorException parse_min(file)
        write(file, "p MDFP 2 0 -1\n")
        @test_throws ErrorException parse_min(file)
        write(file, "p min 2 1\na 1 3 0 1 0\n")
        @test_throws ErrorException parse_min(file)
        write(file, "p min 2 0\nn 1 1\nn 1 -1\n")
        @test_throws ErrorException parse_min(file)
    end
end
