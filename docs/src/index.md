# DiFloT.jl

**DiFloT.jl** is a high-performance Julia package for solving the discrete-time Maximum Dynamic Flow Problem (MDFP) and Quickest Flow Problem (QFP).

It implements:
- The temporally repeated flow algorithm of Ford-Fulkerson.
- Parametric search methods (Binary, Improved Binary, Interpolation, Newton).
- The Lin-Jaillet parametric LP formulation.

## Quick Start

```julia
using DiFloT

# 1. Load a network from a DIMACS-style .min file
net = parse_min("network.min")

# 2. Solve Maximum Dynamic Flow (MDFP)
mdf = solve_max_dynamic_flow(net, 1, 5, 10)
println("Max Flow: ", mdf.value)

# 3. Solve Quickest Flow Problem (QFP)
qf = solve_quickest_flow(net, 1, 5, 50; method="INT-SG", max_iter=200, tol=1e-10)
println("Quickest Time: ", qf.time_T)