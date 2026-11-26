using DiFloT
using Printf

function benchmark_qfp(data_dir=joinpath(@__DIR__, "..", "data"); solver="highs", max_flow_solver="ortools")
    methods = lowercase(solver) == "ortools" ?
        ("BS", "IBS", "INT", "INT-SG") :
        ("BS", "IBS", "INT", "INT-SG", "NEWT", "LP")
    demo = NetworkInstance(
        [1,1,1,2,2,3,3,4], [2,3,4,3,5,4,5,5], [2,2,3,1,1,2,2,4], [1,3,7,1,5,1,3,1]
    )
    
    for method in methods
        solve_quickest_flow(
            demo,1,5,12;
            method=method,solver=solver, max_flow_solver=max_flow_solver,with_paths=false
        )
    end

    for file in sort(filter(f -> endswith(f,".qfp"),readdir(data_dir)))
        net = parse_min(joinpath(data_dir,file))
        s,t,F = net.source,net.sink,Int(net.target_value)
        expected = nothing
        for method in methods
            elapsed = @elapsed r = solve_quickest_flow(
                net,s,t,F;
                method=method,solver=solver,max_flow_solver=max_flow_solver,
                solver_method=method == "LP" ? "dual" : "primal", with_paths=false
            )
            
            @assert r.status == STATUS_OPTIMAL
            
            isnothing(expected) && (expected = r.time_T)
            
            @assert r.time_T == expected
            @printf("%s | %-6s | T*=%.0f | %.6fs | MCCP %d\n", file,method,r.time_T,elapsed,r.mccp_calls)
        end
    end
end

benchmark_qfp()
