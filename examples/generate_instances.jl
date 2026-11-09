using DiFloT

function generate_instances()
    files = ["diflot_netgen_1_10.min", "diflot_netgen_1_100.min"]
    data_dir = joinpath(@__DIR__, "..", "data")
    
    for file in files
        path = joinpath(data_dir, file)
        if !isfile(path)
            println("Skipping $file (not found)")
            continue
        end
        
        net = parse_min(path)
        s, t = net.source, net.sink
        println("Processing $file (Source: $s, Sink: $t)")
        
        # 1. Generate MDFP instances (Factors: 3, 10, 30 * Shortest Path)
        sp = compute_shortest_path(net, s, t)
        for factor in [3, 10, 30]
            net_mdfp = deepcopy(net)
            net_mdfp.problem_type = :MDFP
            net_mdfp.target_value = factor * sp
            out_name = replace(file, ".min" => "_T$(factor)SP.mdfp")
            save_min(net_mdfp, joinpath(data_dir, out_name))
            println("  Generated $out_name (T = $(net_mdfp.target_value))")
        end
        
        # 2. Generate QFP instances (Factors: 2, 5, 10 * Max Flow)
        mf = compute_max_static_flow(net, s, t).value
        for factor in [2, 5, 10]
            net_qfp = deepcopy(net)
            net_qfp.problem_type = :QFP
            net_qfp.target_value = factor * mf
            out_name = replace(file, ".min" => "_F$(factor)MF.qfp")
            save_min(net_qfp, joinpath(data_dir, out_name))
            println("  Generated $out_name (F = $(net_qfp.target_value))")
        end
    end
end

generate_instances()