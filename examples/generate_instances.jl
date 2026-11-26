using DiFloT

# Small illustrative inputs only. The paper uses the separate experiment project.
function generate_instances(data_dir=joinpath(@__DIR__,"..","data"))
    for file in sort(filter(f->endswith(f,".min"),readdir(data_dir)))
        net = parse_min(joinpath(data_dir,file))
        s,t = net.source,net.sink
        sp = compute_shortest_path(net,s,t)
        mf = compute_max_static_flow(net,s,t).value
        if !isfinite(sp) || mf == 0
            println("Skipping disconnected network $file")
            continue
        end
        for (kind,factors,scale,suffix) in ((:MDFP,(3,10,30),Int(sp),"SP"),(:QFP,(2,5,10),mf,"MF"))
            for factor in factors
                out = deepcopy(net)
                out.problem_type = kind
                out.target_value = Float64(Base.Checked.checked_mul(factor,scale))
                name = replace(file,".min"=>"_$(kind)_$(factor)$(suffix).$(lowercase(string(kind)))")
                path = joinpath(data_dir,name)
                if !isfile(path)
                    save_min(out,path)
                    println("Generated $name")
                end
            end
        end
    end
end

generate_instances()
