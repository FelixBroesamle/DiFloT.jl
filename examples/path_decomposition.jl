using DiFloT

# Running example from the manuscript
starts = [1, 1, 1, 2, 2, 3, 3, 4]
ends   = [2, 3, 4, 3, 5, 5, 4, 5]
caps   = [2, 2, 3, 1, 1, 2, 2, 4]
tts    = [1, 3, 7, 1, 5, 3, 1, 1]

s, t = 1, 5
net = NetworkInstance(starts, ends, caps, tts; source=s, sink=t)

# Display a temporally repeated flow
function print_paths(net, paths, T)
    for (i, (P, rate)) in enumerate(paths)
        nodes = [net.arcs[first(P)].tail; [net.arcs[a].head for a in P]]

        transit_time = sum(net.arcs[a].transit_time for a in P)
        value = rate * (T + 1 - transit_time)

        println("Path $i: ", join(nodes, " -> "), " | Rate: ", rate, " | Transit time: ", transit_time, " | Value: ", value)
    end
end

# 1. Maximum Dynamic Flow Problem
T = 6

mdf = solve_max_dynamic_flow(
    net, s, t, T;
    method="FF", solver="HiGHS", with_paths=true
)

println("\n--- MDFP (T = $T) ---")
println("Maximum dynamic flow: ", mdf.value)
print_paths(net, mdf.paths, T)

# 2. Quickest Flow Problem
F = 7

qf = solve_quickest_flow(
    net, s, t, F;
    method="INT-SG", solver="HiGHS", with_paths=true
)

println("\n--- QFP (F = $F) ---")
println("Minimum time-horizon: ", qf.time_T)
print_paths(net, qf.paths, Int(qf.time_T))
