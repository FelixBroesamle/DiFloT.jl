using Test
using DiFloT

@testset "DiFloT.jl" begin
    
    @testset "Types & constants" begin
        @test STATUS_OPTIMAL    == 0
        @test STATUS_INFEASIBLE == 1
        e = Arc(1, 2, 0, 10, 1)
        @test (e.tail, e.head, e.low, e.up, e.transit_time, e.parallel_index) ==
              (1, 2, 0, 10, 1, 0)
    end

    @testset "NetworkInstance (array ctor)" begin
        net = NetworkInstance([1, 1, 2, 3], [2, 3, 4, 4],
                              [10, 5, 10, 10], [1, 2, 1, 1])
        @test net.num_nodes == 4
        @test net.num_arcs  == 4
        @test length(net.arcs) == 4
        @test length(net.arc_indices[(1, 2)]) == 1
    end

    @testset "parse_min / save_min round-trip" begin
        mktempdir() do dir
            fpath = joinpath(dir, "tiny.mdfp")
            open(fpath, "w") do io
                println(io, "c tiny instance")
                println(io, "p MDFP 3 3 7")
                println(io, "n 1 1")
                println(io, "n 3 -1")
                println(io, "a 1 2 0 10 1")
                println(io, "a 2 3 0 10 1")
                println(io, "a 1 3 0 5 3")
            end
            net = parse_min(fpath)
            @test net.problem_type == :MDFP
            @test net.source == 1
            @test net.sink   == 3
            @test net.num_nodes == 3
            @test net.num_arcs  == 3
            @test net.target_value == 7.0

            out = joinpath(dir, "tiny_out")
            save_min(net, out)
            @test isfile(out * ".mdfp")
            net2 = parse_min(out * ".mdfp")
            @test net2.num_arcs  == net.num_arcs
            @test net2.num_nodes == net.num_nodes
            @test net2.target_value == net.target_value
        end
    end

    @testset "compute_shortest_path" begin
        net = NetworkInstance([1, 1, 2], [2, 3, 3],
                              [5, 5, 5], [1, 5, 1])
        @test compute_shortest_path(net, 1, 3) == 2.0   # 1→2→3
        @test compute_shortest_path(net, 1, 1) == 0.0
    end

    # Diamond: 1→2→4 and 1→3→4, each cap 10, tt 1.
    # MDFP(T) = 20(T+1) - 40 = 20T - 20 for T ≥ 2.
    # QFP(F=20)  = 2.
    let starts = [1, 1, 2, 3], ends = [2, 3, 4, 4],
        caps   = [10, 10, 10, 10], tts = [1, 1, 1, 1]

        net = NetworkInstance(starts, ends, caps, tts)

        @testset "MDFP method=FF solver=HiGHS" begin
            res = solve_max_dynamic_flow(net, 1, 4, 5;
                                         method="FF", solver="HiGHS")
            @test res.status == STATUS_OPTIMAL
            @test res.value  ≈ 80.0 atol = 1e-6
            @test !isempty(res.paths)
        end

        @testset "MDFP method=TE solver=HiGHS" begin
            res = solve_max_dynamic_flow(net, 1, 4, 5;
                                         method="TE", solver="HiGHS")
            @test res.status == STATUS_OPTIMAL
            @test res.value  ≈ 80.0 atol = 1e-6
        end

        # QFP methods that go through the parametric oracle
        for m in ("BS", "IBS", "INT", "INT-SG", "NEWT", "LP")
            @testset "QFP method=$m" begin
                res = solve_quickest_flow(net, 1, 4, 20;
                                          method=m, solver="HiGHS")
                @test res.status == STATUS_OPTIMAL
                @test res.time_T == 2.0
            end
        end
    end

    let starts = [1, 1, 2, 3], ends = [2, 3, 4, 4],
        caps   = [10, 10, 10, 10], tts = [1, 1, 1, 1]

        net = NetworkInstance(starts, ends, caps, tts)
        @testset "Newton correctness" begin
            # --- (a) F = 20: early return at T_l = T* = 2 -------------------------
            #     v(2) = 20, so T_l = 2 and v(T_l) >= F, early exit.
            res = solve_quickest_flow(net, 1, 4, 20; method="NEWT", solver="HiGHS")
            @test res.status == STATUS_OPTIMAL
            @test res.time_T == 2.0

            # --- (b) F = 21: Newton step, then two-step verification -------------
            #     T_l = 2, v(2) = 20 < 21.
            #     θ_1 = 2 + 1/20 = 2.05, v(2.05) = 21, |v - F| = 0 < 1, exit loop.
            #     ⌊2.05⌋ = 2, v(2) = 20 < 21, so return 3.
            res = solve_quickest_flow(net, 1, 4, 21; method="NEWT", solver="HiGHS")
            @test res.status == STATUS_OPTIMAL
            @test res.time_T == 3.0

            # --- (c) F = 39: large Newton step, exit near T* = 3 -----------------
            #     T_l = 2, v(2) = 20 < 39.
            #     θ_1 = 2 + 19/20 = 2.95, v(2.95) = 39, exit.
            #     ⌊2.95⌋ = 2, v(2) < 39, return 3.
            res = solve_quickest_flow(net, 1, 4, 39; method="NEWT", solver="HiGHS")
            @test res.status == STATUS_OPTIMAL
            @test res.time_T == 3.0

            # --- (d) F = 40: exits on integer θ with |v - F| = 0 ------------------
            #     T_l = 2, θ_1 = 3.0 exactly, v(3) = 40, exit.
            #     ⌊3⌋ = 3, v(3) = 40 >= 40, return 3.
            res = solve_quickest_flow(net, 1, 4, 40; method="NEWT", solver="HiGHS")
            @test res.status == STATUS_OPTIMAL
            @test res.time_T == 3.0

            # --- (e) Consistency with BS on all F in a small sweep ---------------
            for F in 20:1:80
                r_bs   = solve_quickest_flow(net, 1, 4, F; method="BS",   solver="HiGHS")
                r_newt = solve_quickest_flow(net, 1, 4, F; method="NEWT", solver="HiGHS")
                @test r_newt.time_T == r_bs.time_T
                @test r_newt.status == STATUS_OPTIMAL
            end

            # --- (f) Infeasible instance: g_max = 0 -------------------------------
            # Arc 1→2 has zero capacity; no flow can leave s.
            net_inf = NetworkInstance([1, 2], [2, 4], [0, 10], [1, 1])
            res = solve_quickest_flow(net_inf, 1, 4, 5; method="NEWT", solver="HiGHS")
            @test res.status == STATUS_INFEASIBLE
        end
    end
end

