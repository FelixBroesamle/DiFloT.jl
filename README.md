# DiFloT.jl

[![CI](https://github.com/FelixBroesamle/DiFloT.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/FelixBroesamle/DiFloT.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/FelixBroesamle/DiFloT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/FelixBroesamle/DiFloT.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Julia Version](https://img.shields.io/badge/julia-1.10%2B-blue)](https://julialang.org/)

**A high-performance Julia solver for the Maximum Dynamic Flow Problem (MDFP) and Quickest Flow Problem (QFP) in Discrete Time.**

**DiFloT.jl** (**Di**screte **Flo**ws over **T**ime) is designed for researchers in Operations Research, Network Optimization, and Mathematical Programming and for researchers and practitioners in areas of applications, such as Evacuation Planning, Supply Chain Optimizatoin, and Data Transmission. It provides highly optimized, state-of-the-art algorithms for solving dynamic flow problems over discrete time horizons.

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

### Standard installation
```julia
using Pkg
Pkg.add("https://github.com/FelixBroesamle/DiFloT.jl")
```

### Building from Source
```bash
git clone https://github.com/FelixBroesamle/DiFloT.jl.git
cd DiFloT.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Verifying the Installation
```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Optional Gurobi Support

By default, `DiFloT.jl` uses the open-source HiGHS solver. If you have a Gurobi license and wish to use it, simply add and load Gurobi in your environment:

```julia
using Pkg; Pkg.add("Gurobi")
using DiFloT, Gurobi
# Now you can pass solver="Gurobi"
```

### Running the Example Workflows

The repository includes demo scripts and datasets in the `data/` and `examples/` directories.

```bash
# Run the quickstart tutorial
julia --project=. examples/quickstart.jl

# Generate .mdfp and .qfp instances from base .min files
julia --project=. examples/generate_instances.jl

# Run (small / illustrative) benchmarks for MDFP and QFP
julia --project=. examples/benchmark_mdfp.jl
julia --project=. examples/benchmark_qfp.jl
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

# 2. Solve Maximum Dynamic Flow (MDFP)
s = 1
t = 4
time_horizon = 10
mdf = solve_max_dynamic_flow(net, s, t, time_horizon; method="FF")
println("Max Flow: ", mdf.value)

# 3. Solve Quickest Flow Problem (QFP)
f_value = 20
qf = solve_quickest_flow(net, s, t, f_value; method="IBS")
println("Quickest Time: ", qf.time_T)
```

## File Formats

`DiFloT.jl` uses a DIMACS-style format, separated by file extension for the two problem classes:
- `.min`: Standard DIMACS format.
- `.mdfp`: Extension to Maximum Dynamic Flow instances. Includes the time horizon `T` in the `p` line.
- `.qfp`: Extension to Quickest Flow instances. Includes the flow value `F` in the `p` line.

The `p` line defines the problem type, dimensions and problem parameter (`T`, the time-horizon for the MDFP, and `F`, the flow value for the QFP). The `n` lines define the source (`1`) and sink (`-1`). The `a` lines represent `(tail, head, lower_bound, capacity, transit_time)`.

```text
c This is an MDFP instance
p MDFP 3 3 7
n 1 1
n 3 -1
a 1 2 0 10 1
a 2 3 0 10 1
a 1 3 0 5 3
```
```text
c This is a QFP instance
p MDFP 3 3 42
n 1 1
n 3 -1
a 1 2 0 10 1
a 2 3 0 10 1
a 1 3 0 5 3
```

## End-to-End Optimization Workflows

The file-based solvers allow you to load instances directly from disk. `DiFloT.jl` will automatically parse the parameter `T` or `F` from `.mdfp` / `.qfp` files.

### Workflow 1: From File (Auto-detect parameters)

When loading a `.mdfp` or `.qfp` file, `DiFloT.jl` automatically extracts the source, sink, and target values (`T` or `F`) directly from the file headers.

```julia
using DiFloT

# Loads the problem and automatically extracts T from the file
res = solve_max_dynamic_flow_from_file("data/network.mdfp", s=1, t=3)

# Loads the problem and automatically extracts F from the file
res_qf = solve_quickest_flow_from_file("data/network.qfp", s=1, t=3)
```

### Workflow 2: From File (Override parameters)

If you load a `.min` file (and interpret the costs as transit times) or want to override the values in the file, simply pass `T` or `F` explicitly:

```julia
# Overrides the time horizon in the .mdfp file with T=50
res = solve_max_dynamic_flow_from_file("data/network.mdfp", 1, 3, 50)

# Solves a .min file by providing T manually
res_min = solve_max_dynamic_flow_from_file("data/base.min", 1, 3, 20)
```

### Workflow 3: In-Memory Generation

For integration into larger pipelines, you can build a `NetworkInstance` directly from arrays and pass it to the solvers.

```julia
using DiFloT

# Define a network directly via arrays (Tail, Head, Capacity, Transit Time)
starts = [1, 1, 2, 3]
ends   = [2, 3, 4, 4]
caps   = [10, 5, 10, 10]
tts    = [1, 2, 1, 1]

# 1. Build the NetworkInstance
net = NetworkInstance(starts, ends, caps, tts)

# 2. Solve MDFP (T=10)
mdf = solve_max_dynamic_flow(net, 1, 4, 10; method="FF")

# 3. Solve QFP (F=20)
qf = solve_quickest_flow(net, 1, 4, 20; method="IBS")
```

## Citing

If you use `DiFloT.jl` in your research, please cite the following paper:

```text
To appear.
```

## Resources

- Documentation: (To appear). [https://FelixBroesamle.github.io/DiFloT.jl/](https://FelixBroesamle.github.io/DiFloT.jl/).
- Package: [https://github.com/FelixBroesamle/DiFloT.jl](https://github.com/FelixBroesamle/DiFloT.jl).

## License

Distributed under the MIT License. See `LICENSE` for more information.

