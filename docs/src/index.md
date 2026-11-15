# DiFloT.jl

**A high-performance Julia solver for the Maximum Dynamic Flow Problem (MDFP) and Quickest Flow Problem (QFP) in Discrete Time.**

**DiFloT.jl** (**Di**screte **Flo**ws over **T**ime) is designed for researchers in Operations Research, Network Optimization, and Mathematical Programming, as well as for practitioners in areas such as Evacuation Planning, Supply Chain Optimization, and Data Transmission. It provides highly optimized, state-of-the-art algorithms for solving dynamic flow problems over discrete time horizons.

## Key Features

- **Maximum Dynamic Flow Problem (MDFP)**:
    - The Ford-Fulkerson Algorithm (`"FF"`).
    - Solving the Problem on the time-expanded network (`"TE"`).
- **Quickest Flow Problem (QFP)**:
    - Standard binary search (`"BS"`).
    - Improved binary search (`"IBS"`).
    - Interpolation search (`"INT"`) and a safeguarded Interpolation search (`"INT-SG"`).
    - Newton search (`"NEWT"`).
    - Lin-Jaillet LP formulation (`"LP"`).
- **DIMACS Compatible**: Parses `.min` networks and introduces `.mdfp` and `.qfp` formats to standardize dynamic flow instances.


## Installation

### Standard Installation
`DiFloT.jl` can be installed directly from the Julia General Registry:

```julia
using Pkg
Pkg.add(url="https://github.com/FelixBroesamle/DiFloT.jl")
```

### Optional Gurobi Support

By default, `DiFloT.jl` uses the open-source HiGHS solver. If you have a Gurobi license and wish to use it, simply add and load Gurobi in your environment:

```julia
using Pkg; Pkg.add("Gurobi")
using DiFloT, Gurobi
# Now you can pass solver="Gurobi"
```

## Quick Start

```julia
using DiFloT

# 1. Define a network directly via arrays
starts = [1, 1, 2, 3]
ends   = [2, 3, 4, 4]
caps   = [10, 5, 10, 10]
tts    = [1, 2, 1, 1]
net = NetworkInstance(starts, ends, caps, tts)

# 2. Solve Maximum Dynamic Flow (MDFP) with time horizon T=10
mdf = solve_max_dynamic_flow(net, 1, 4, 10; method="FF")
println("Max Flow (MDFP) at T=10: ", mdf.value)

# 3. Solve Quickest Flow Problem (QFP) for target flow F=20
qf = solve_quickest_flow(net, 1, 4, 20; method="IBS")
println("Quickest Time (QFP) for F=20: ", qf.time_T)
```




