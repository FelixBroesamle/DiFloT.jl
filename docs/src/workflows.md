# End-to-End Workflows

The repository includes demo scripts and datasets in the `data/` and `examples/` directories. This page highlights the primary workflows. To run them, ensure you have the repository cloned locally.

## 1. In-Memory Generation (Array Input)

For direct integration into larger Julia pipelines, you can build a `NetworkInstance` directly from arrays. This bypasses file I/O entirely.

```julia
using DiFloT

# 1. Define a network via arrays (Tail, Head, Capacity, Transit Time)
starts = [1, 1, 2, 3]
ends   = [2, 3, 4, 4]
caps   = [10, 5, 10, 10]
tts    = [1, 2, 1, 1]
net = NetworkInstance(starts, ends, caps, tts)

# 2. Solve Maximum Dynamic Flow (MDFP) for a time horizon T=10
mdf = solve_max_dynamic_flow(net, 1, 4, 10; method="FF")
println("Max Flow (MDFP) at T=10: ", mdf.value)

# 3. Solve Quickest Flow Problem (QFP) for a demand F=20
qf = solve_quickest_flow(net, 1, 4, 20; method="IBS")
println("Quickest Time (QFP) for F=20: ", qf.time_T)
```

## 2. File-Based Parsing (Auto-detect parameters)

When loading a `.mdfp` or `.qfp` file, DiFloT.jl automatically extracts the source, sink, and target values (`T` or `F`) directly from the file headers.

```julia
using DiFloT

# 1. Solve MDFP by loading a .mdfp file
# Automatically reads s, t, and T from the file
res_mdfp = solve_max_dynamic_flow_from_file("data/network.mdfp", method="FF")
println("Max Flow: ", res_mdfp.value)

# 2. Solve QFP by loading a .qfp file
# Automatically reads s, t, and F from the file
res_qfp = solve_quickest_flow_from_file("data/network.qfp", method="INT-SG")
println("Quickest Time: ", res_qfp.time_T)
```

## 3. File-Parsing (Override parameters)

If you load a standard `.min` file (which lacks the input parameters `T` or `F`), or if you want to overwrite the values embedded in a `.mdfp` or `.qfp` file, you can pass `T` or `F` explicitly.

```julia
using DiFloT

# 1. Load a base .min file. s and t are auto-detected, but T must be provided.
res_min_mdfp = solve_max_dynamic_flow_from_file("data/base_network.min", T=50)
println("Max Flow (T=50): ", res_min_mdfp.value)

# 2. Load a .qfp file but override the demand to F=100
res_min_qfp = solve_quickest_flow_from_file("data/network.qfp", F=100)
println("Quickest Time (F=100): ", res_min_qfp.time_T)
```
