using DiFloT

println("--- DiFloT.jl Quickstart ---")

# 1. Define a network directly via arrays
starts = [1, 1, 2, 3]
ends   = [2, 3, 4, 4]
caps   = [10, 5, 10, 10]
tts    = [1, 2, 1, 1]
net = NetworkInstance(starts, ends, caps, tts)

# 2. Solve Maximum Dynamic Flow (MDFP) with time horizon T=10
s = 1
t = 4
time_horizon = 10
mdf = solve_max_dynamic_flow(net, s, t, time_horizon; method="FF")
println("Max Flow (MDFP) at T=10: ", mdf.value)

# 3. Solve Quickest Flow Problem (QFP) for target flow F=20
f_value = 20
qf = solve_quickest_flow(net, s, t, f_value; method="IBS")
println("Quickest Time (QFP) for F=20: ", qf.time_T)