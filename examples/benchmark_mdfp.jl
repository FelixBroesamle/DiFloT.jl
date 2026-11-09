using DiFloT
using Printf

function benchmark_mdfp()
    data_dir = joinpath(@__DIR__, "..", "data")
    files = filter(f -> endswith(f, ".mdfp"), readdir(data_dir))
    
    if isempty(files)
        println("No .mdfp files found. Run examples/generate_instances.jl first!")
        return
    end
    
    for file in files
        path = joinpath(data_dir, file)
        net = parse_min(path)
        s, t_node, T = net.source, net.sink, Int(net.target_value)
        
        @printf("\n--- Benchmarking %s (T=%d) ---\n", file, T)
        
        # 1. Ford-Fulkerson (MCCP) with HiGHS
        t_ff = @elapsed res_ff = solve_max_dynamic_flow(net, s, t_node, T; method="FF", solver="HiGHS")
        @printf("FF (HiGHS MCCP)        : %.3fs | Value: %.1f\n", t_ff, res_ff.value)
        
        # 2. Time-Expanded with pure Julia Sparse Dinic
        t_te = @elapsed res_te = solve_max_dynamic_flow(net, s, t_node, T; method="TE", solver="GraphsFlows")
        @printf("TE (Sparse Dinic)      : %.3fs | Value: %.1f\n", t_te, res_te.value)
        
        # 3. Time-Expanded with HiGHS LP
        #t_lp = @elapsed res_lp = solve_max_dynamic_flow(net, s, t_node, T; method="TE", solver="HiGHS")
        #@printf("TE (HiGHS LP)          : %.3fs | Value: %.1f\n", t_lp, res_lp.value)
    end
end

benchmark_mdfp()