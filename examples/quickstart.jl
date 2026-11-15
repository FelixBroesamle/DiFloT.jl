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

# 4. Example from manuscript
starts_exp = [1, 1, 1, 2, 2, 3, 3, 4]
ends_exp = [2, 3, 4, 3, 5, 4, 5, 5]
caps_exp = [2, 2, 3, 1, 1, 2, 2, 4]
tts_exp = [1, 3, 7, 1, 5, 1, 3, 1]

net_exp = NetworkInstance(starts_exp, ends_exp, caps_exp, tts_exp)

s_exp = 1
t_exp = 5
time_horizon_exp_1 = 9
time_horizon_exp_2 = 15

mdf_exp_1 = solve_max_dynamic_flow(net_exp, s_exp, t_exp, time_horizon_exp_1; method="FF")
mdf_exp_2 = solve_max_dynamic_flow(net_exp, s_exp, t_exp, time_horizon_exp_2; method="FF")
println("Max Flow (MDFP) at T=9: ", mdf_exp_1.value)
println("Max Flow (MDFP) at T=15: ", mdf_exp_2.value)


net_exp = NetworkInstance(starts_exp, ends_exp, caps_exp, tts_exp)
s_exp, t_exp = 1, 5

println("T  | v(T)")
println("---+-----")

for T in 1:16
    res = solve_max_dynamic_flow(net_exp, s_exp, t_exp, T; method="FF")
    v_T = Int(res.value)
    println(lpad(T, 2), " | ", lpad(v_T, 4))
end
