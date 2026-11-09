using DiFloT
using Printf

function benchmark_qfp()
    data_dir = joinpath(@__DIR__, "..", "data")
    files = filter(f -> endswith(f, ".qfp"), readdir(data_dir))
    methods = ["BS", "IBS", "INT", "INT-SG", "NEWT", "LP"]
    
    for file in files
        path = joinpath(data_dir, file)
        net = parse_min(path)
        
        # RENAMED `t` to `t_node`
        s, t_node, F = net.source, net.sink, net.target_value
        
        @printf("\n--- Benchmarking %s (F=%.1f) ---\n", file, F)
        
        for m in methods
            # Use `t_node` here so the elapsed time doesn't overwrite it
            elapsed = @elapsed res = solve_quickest_flow(net, s, t_node, F; method=m, solver="HiGHS")
            @printf("%-6s : %.3fs | Time T*: %.1f\n", m, elapsed, res.time_T)
        end
    end
end

benchmark_qfp()