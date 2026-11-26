@testset "OR-Tools and HiGHS integration" begin
    # Deliberately solve with HiGHS BEFORE the worker imports OR-Tools.
    net = NetworkInstance([1,1],[2,2],[3,4],[0,1])
    @test solve_max_dynamic_flow(net,1,2,2;solver="HiGHS",solver_method="primal").value == 17
    version = start_ortools!() # No catch or skip: a requested backend must work.
    @test !isempty(version)
    worker = DiFloT._ORTOOLS_IO[]
    @test worker !== nothing
    @test start_ortools!() == version
    @test DiFloT._ORTOOLS_IO[] === worker
    @info "Testing real OR-Tools worker" version Sys.ARCH

    try
        @testset "Both real OR-Tools graph solvers" begin
            # Force the raw IPC path; public methods can legitimately return early.
            raw = DiFloT._ortools_request(1,[1,1],[2,2],[3,4],2,1,2)
            @test raw.value == 7
            @test isempty(raw.flows)
            
            # Zero supplies, negative return-arc cost: minimum-cost circulation.
            raw = DiFloT._ortools_request(
                2,[1,1,2],[2,2,1],[3,4,7],2,1,2;
                costs=[0,1,-3]
            )
            @test raw.value == 17
            @test raw.flows == [3,4,7]
            
            # Disconnected components, with nonzero terminal capacity bounds.
            # This cannot be answered by a trivial zero-rate precheck.
            disconnected = NetworkInstance([1,3],[2,4],[5,5],[1,1];source=1,sink=4)
            @test DiFloT.c_rate_bound(disconnected,1,4) == 5
            mf = compute_max_static_flow(disconnected,1,4;solver="ORTools")
            @test mf.status == STATUS_OPTIMAL && mf.value == 0
        end

        @testset "MDFP: three independent configurations" begin
            # Positive/zero horizons, parallel arcs and a zero-time arc.
            for T in 0:4
                expected = 7*T+3
                hi = solve_max_dynamic_flow(net,1,2,T;solver="HiGHS",solver_method="dual")
                ff = solve_max_dynamic_flow(net,1,2,T;solver="ORTools")
                te = solve_max_dynamic_flow(net,1,2,T;method="TE",max_flow_solver="ORTools")
                @test hi.status == ff.status == te.status == STATUS_OPTIMAL
                @test hi.value == ff.value == te.value == expected
                @test repeated_value(net,T,ff.paths) == expected
                check_paths(net,1,2,T,ff.paths)
            end
            
            # A nontrivial network with shared arcs and several breakpoints.
            example = parse_min(joinpath(@__DIR__,"..","data","example.min"))
            for (T,expected) in ((3,0),(4,1),(5,3),(6,7),(7,11),(8,17),(15,66))
                hi = solve_max_dynamic_flow(example,1,5,T;solver="HiGHS",solver_method="dual",with_paths=false)
                ff = solve_max_dynamic_flow(example,1,5,T;solver="ORTools",with_paths=false)
                te = solve_max_dynamic_flow(example,1,5,T;method="TE",max_flow_solver="ORTools",with_paths=false)
                @test hi.status == ff.status == te.status == STATUS_OPTIMAL
                @test hi.value == ff.value == te.value == expected
                @test isempty(hi.paths) && isempty(ff.paths)
            end
        end

        @testset "Directed and parallel arcs" begin
            directed_net = NetworkInstance(
                [1,2,1,2], [2,1,2,3], [3,99,5,8], [0,7,2,1]; 
                source=1, sink=3
            )
            @test compute_max_static_flow(directed_net, 1, 3; solver="ORTools").value == 8
            @test compute_max_static_flow(directed_net, 3, 1; solver="ORTools").value == 0
            
            for T in 1:5
                expected = 3*T + 5*max(0,T-2)
                ff = solve_max_dynamic_flow(directed_net, 1, 3, T; method="FF", solver="ORTools")
                te = solve_max_dynamic_flow(directed_net, 1, 3, T; method="TE", max_flow_solver="ORTools")
                @test ff.value == te.value == expected
                check_paths(directed_net, 1, 3, T, ff.paths)
            end
        end

        @testset "QFP: mixed-backend configuration" begin
            example = parse_min(joinpath(@__DIR__,"..","data","example.min"))
            
            # Known optima, independently checked via OR-Tools time expansion.
            for (F,T) in ((1,4),(7,6),(12,8),(24,9),(25,10))
                v = solve_max_dynamic_flow(example,1,5,T;method="TE",max_flow_solver="ORTools").value
                prev = solve_max_dynamic_flow(example,1,5,T-1;method="TE",max_flow_solver="ORTools").value
                @test prev < F <= v
            
                for method in ("BS","IBS","INT","INT-SG","NEWT")
            
                    r = solve_quickest_flow(
                        example,1,5,F;
                        method=method,solver="HiGHS", solver_method="primal",max_flow_solver="ORTools",
                        with_paths=false,verify=false
                    )
                    @test r.status == STATUS_OPTIMAL && r.time_T == T
                    @test r.mccp_calls >= 1 && r.lp_calls == 0
                end
                
                candidate = solve_lp_candidate(example,1,5,F;solver="HiGHS",solver_method="dual")
                @test candidate.status == STATUS_OPTIMAL && candidate.lp_calls == 1
                @test abs(candidate.time_T-T) <= 1 # LP-only is not certified exact.
                
                checked = check_lp_candidate(
                    example,1,5,F,candidate;
                    solver="HiGHS", solver_method="primal",with_paths=false,verify=false
                )
                @test checked.status == STATUS_OPTIMAL && checked.time_T == T
                @test 1 <= checked.mccp_calls <= 2 && checked.lp_calls == 1
            end
        end
        
        @testset "Integer QFP oracles with OR-Tools" begin
            example = parse_min(joinpath(@__DIR__, "..", "data", "example.min"))
            for (F, T) in ((1,4), (7,6), (12,8), (24,9), (25,10))
                for method in ("BS", "IBS", "INT", "INT-SG")
                    r = solve_quickest_flow(
                        example, 1, 5, F;
                        method=method, solver="ORTools", max_flow_solver="ORTools"
                    )
                    @test r.status == STATUS_OPTIMAL && r.time_T == T
                    @test r.mccp_calls >= 1
                end
            end
            for method in ("NEWT", "LP")
                @test_throws ArgumentError solve_quickest_flow(
                    example, 1, 5, 12;
                    method=method, solver="ORTools"
                )
            end
            session = DiFloT.MCCPSession(net, 1, 2; solver="ORTools")
            @test DiFloT.solve_mccp(session, 2).v == 17
            @test DiFloT.solve_mccp(session, 2).v == 17
            @test session.calls == 1
            @test_throws ArgumentError DiFloT.solve_mccp(session, 2.5)
            @test DiFloT.solve_mccp(session, 0).v == 3
            @test session.calls == 2
        end
        @test DiFloT._ORTOOLS_IO[] === worker # Repeated requests reused the same process.
        
        # HiGHS is still usable after all OR-Tools operations.
        @test solve_max_dynamic_flow(net,1,2,2;solver="HiGHS",solver_method="primal").value == 17
    finally
        stop_ortools!()
    end
    @test DiFloT._ORTOOLS_IO[] === nothing
    stop_ortools!() # Idempotent shutdown.

    # Starting from a stopped worker proves the QFP actually requests OR-Tools.
    try
        r = solve_quickest_flow(
            net,1,2,8;
            method="BS",solver="HiGHS", solver_method="primal",max_flow_solver="ORTools"
        )
        @test r.status == STATUS_OPTIMAL && r.time_T == 1
        @test DiFloT._ORTOOLS_IO[] !== nothing
        @test DiFloT._ORTOOLS_IO[] !== worker
        @test start_ortools!() == version
    finally
        stop_ortools!()
    end
end
