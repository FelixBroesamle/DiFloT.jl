using Test
using DiFloT

# Optional licensed test: julia --project=. test/optional_backends.jl Gurobi
if "Gurobi" in ARGS
    using Gurobi
    
    @testset "Optional Gurobi" begin
        net = NetworkInstance([1,1],[2,2],[3,4],[0,1])
    
        for sm in ("primal","dual")
            @test solve_max_dynamic_flow(net,1,2,2;solver="Gurobi",solver_method=sm).value == 17
    
            for method in ("BS","IBS","INT","INT-SG","NEWT","LP")
                @test solve_quickest_flow(net,1,2,8;method=method,solver="Gurobi",solver_method=sm,verify=true).time_T == 1
            end
        end
    end
end
