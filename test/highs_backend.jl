@testset "HiGHS simplex strategies and reoptimization" begin
    # v(theta) = 3*(theta+1) + 4*theta + 2*max(0,theta-3).
    net = NetworkInstance([1,1,1],[2,2,2],[3,4,2],[0,1,4])
    
    for (method,strategy) in (("primal",4),("dual",1))
        @testset "$method" begin
            session = DiFloT.MCCPSession(net,1,2;solver="HiGHS",solver_method=method)
            model = session.model
    
            for theta in (0,2,0.5,4,1,2)
                r = DiFloT.solve_mccp(session,theta)
                @test r.status == STATUS_OPTIMAL
                @test r.v ≈ 7*theta+3+2*max(0,theta-3) atol=1e-8
                @test session.model === model
    
                if theta > 0 && theta != 3
                    @test r.g == (theta < 3 ? 7 : 9)
                end
    
                # Compare to a new session with no prior solves.
                cold = DiFloT.MCCPSession(net,1,2;solver="HiGHS",solver_method=method)
                @test DiFloT.solve_mccp(cold,theta).v ≈ r.v atol=1e-8
            end
    
            @test JuMP.get_optimizer_attribute(model,"solver") == "simplex"
            @test JuMP.get_optimizer_attribute(model,"simplex_strategy") == strategy
            @test JuMP.get_optimizer_attribute(model,"threads") == 1
            @test session.calls == 6
    
            old = session.last_result
            @test DiFloT.solve_mccp(session,2) === old
            @test session.calls == 6
            @test DiFloT.solve_mccp(session,2;force=true).v == 17
            @test session.calls == 7
        end
    end
end

@testset "MCCP cache and interpolation" begin
    net = NetworkInstance([1],[2],[3],[2])
    session = DiFloT.MCCPSession(net,1,2;solver_method="dual")
    r = DiFloT.solve_mccp(session,2)
    
    @test r.v == 3
    @test DiFloT.solve_mccp(session,2) === r
    @test session.calls == 1
    @test DiFloT.solve_mccp(session,2.5).v ≈ 4.5
    @test session.calls == 2
    @test DiFloT.solve_mccp(session,2).v == 3
    @test session.calls == 3
    @test DiFloT._illinois(x->x^2-0.25,-0.25,0.75,100,1e-10) ≈ 0.5 atol=1e-8
    
    # Interpolation example, including a steep exponent b=24.
    r_l = DiFloT.MCCPResult(0,1.0,1.0,Float64[])
    r_u = DiFloT.MCCPResult(0,66.0,7.0,Float64[])
    
    @test DiFloT._interp_candidate(4,15,r_l,r_u,12,100,1e-9) == 7
    r_l = DiFloT.MCCPResult(0,11.0,4.0,Float64[])
    
    @test DiFloT._interp_candidate(7,15,r_l,r_u,12,100,1e-9) == 8
    # Root exactly at upper endpoint; projection still guarantees progress.
    r_l = DiFloT.MCCPResult(0,1.0,1.0,Float64[])
    r_u = DiFloT.MCCPResult(0,11.0,1.0,Float64[])
    
    @test DiFloT._interp_candidate(1,11,r_l,r_u,11,100,1e-9) == 10
end

