using DiFloT

# 1. In-memory network.
net = NetworkInstance([1,1,2,3],[2,3,4,4],[10,5,10,10],[1,2,1,1];source=1,sink=4)
mdf = solve_max_dynamic_flow(net,1,4,10;method="FF")
qf = solve_quickest_flow(net,1,4,20;method="IBS")
@assert mdf.value == 130 && qf.time_T == 3
println("MDFP value: ",mdf.value," | QFP horizon: ",qf.time_T)

# 2. Running example from the manuscript.
net = parse_min(joinpath(@__DIR__,"..","data","example.min"))
s,t = net.source,net.sink
println("T  | v(T)")
for T in 0:16
    r = solve_max_dynamic_flow(net,s,t,T;method="FF",with_paths=false)
    println(lpad(T,2)," | ",lpad(Int(r.value),4))
end
for method in ("BS","IBS","INT","INT-SG","NEWT","LP")
    r = solve_quickest_flow(net,s,t,12;method=method)
    @assert r.time_T == 8
    println(method,": T* = ",Int(r.time_T),", MCCP calls = ",r.mccp_calls)
end
