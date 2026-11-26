using DiFloT
using Printf

# Illustrative timings with default open-source Julia backends.
function benchmark_mdfp(data_dir=joinpath(@__DIR__,"..","data"))
    demo = parse_min(joinpath(@__DIR__,"..","data","example.min"))
    for method in ("FF","TE")
        solve_max_dynamic_flow(demo,1,5,5;method=method,with_paths=false)
    end
    for file in sort(filter(f->endswith(f,".mdfp"),readdir(data_dir)))
        net = parse_min(joinpath(data_dir,file))
        s,t,T = net.source,net.sink,Int(net.target_value)
        elapsed_ff = @elapsed ff = solve_max_dynamic_flow(net,s,t,T;method="FF",solver_method="dual",with_paths=false)
        elapsed_te = @elapsed te = solve_max_dynamic_flow(net,s,t,T;method="TE",solver="ortools",with_paths=false)
        @assert ff.value == te.value
        @printf("%s | T=%d | value=%.0f | FF %.6fs | TE %.6fs\n",file,T,ff.value,elapsed_ff,elapsed_te)
    end
end

benchmark_mdfp()
