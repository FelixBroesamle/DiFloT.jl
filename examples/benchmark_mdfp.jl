using DiFloT
using Printf

function benchmark_mdfp()
    data_dir = joinpath(@__DIR__, "..", "data")
    files = filter(f -> endswith(f, ".mdfp"), readdir(data_dir))
    
    for file in files
        path = joinpath(data_dir, file)
        net = parse_min(path)
        s, t, T = net.source, net.sink, Int(net.target_value)
        
        @printf("\n--- Benchmarking %s (T=%d) ---\n", file, T)
        
        # FF with HiGHS (MCCP)
        t_ff_highs = @elapsed res_ff = solve_max_dynamic_flow(net, s, t, T; method="FF", solver="HiGHS")
        @printf("FF (HiGHS)      : %.3fs | Value: %.1f\n", t_ff_highs, res_ff.value)
        
        # TE with pure Julia GraphsFlows
        #t_te_graph = @elapsed res_te = solve_max_dynamic_flow(net, s, t, T; method="TE")
        #@printf("TE (GraphsFlows): %.3fs | Value: %.1f\n", t_te_graph, res_te.value)
    end
end

benchmark_mdfp()