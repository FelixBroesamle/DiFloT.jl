# DiFloT.jl

[![CI](https://github.com/FelixBroesamle/DiFloT.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/FelixBroesamle/DiFloT.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/FelixBroesamle/DiFloT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/FelixBroesamle/DiFloT.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Julia Version](https://img.shields.io/badge/julia-1.10%2B-blue)](https://julialang.org/)

**A high-performance Julia solver for the Maximum Dynamic Flow Problem (MDFP) and Quickest Flow Problem (QFP) in Discrete Time.**

**DiFloT.jl** (**Di**screte **Flo**ws over **T**ime) is designed for researchers in Operations Research, Network Optimization, and Mathematical Programming and for researchers and practitioners in areas of applications, such as Evacuation Planning, Supply Chain Optimization, and Data Transmission. It provides highly optimized, state-of-the-art algorithms for solving dynamic flow problems over discrete time horizons.

## Key Features

- **Maximum Dynamic Flow Problem (MDFP)**:
    - The Ford-Fulkerson Algorithm (`"FF"`).
    - Solving the Problem on the time-expanded network (`"TE"`).
- **Quickest Flow Problem (QFP)**:
    - Standard binary search (`"BS"`).
    - Improved binary search (`"IBS"`).
    - Interpolation search (`"INT"`).
    - Safeguarded Interpolation search (`"INT-SG"`).
    - Newton search (`"NEWT"`).
    - Lin-Jaillet LP formulation (`"LP"`).
- **Solver backends:** HiGHS and GraphsFlows by default,
  with optional OR-Tools and Gurobi support.
- **Network representation:** Directed networks with
  parallel arcs and nonnegative integer transit times.
- **DIMACS Compatible**: Parses `.min` networks and introduces `.mdfp` and `.qfp` formats to standardize dynamic flow instances.

## Installation

### Standard installation
Install DiFloT.jl directly from GitHub:

```julia-repl
using Pkg
Pkg.add(url="https://github.com/FelixBroesamle/DiFloT.jl")
```

### Building from Source

```bash
git clone https://github.com/FelixBroesamle/DiFloT.jl.git
cd DiFloT.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Verifying the Installation

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The default installation uses HiGHS and GraphsFlows.
No Python installation or commercial solver license is required.


### Optional OR-Tools Support

OR-Tools can be used to solve maximum-flow problems
and minimum-cost circulation problems with integer costs.

To avoid native-library conflicts between OR-Tools and HiGHS,
DiFloT.jl runs OR-Tools in a separate Python process.

The repository provides a `CondaPkg.toml` file specifying
the Python environment and its dependencies.

From the cloned repository, activate the project and prepare
the Python environment:

```julia-repl
using Pkg
Pkg.activate(".")
Pkg.add("CondaPkg")

using CondaPkg
CondaPkg.resolve()

ENV["DIFLOT_PYTHON"] = CondaPkg.which("python")

using DiFloT

CondaPkg.withenv() do
    start_ortools!()
end
```

The worker is reused across successive optimization calls
and can be stopped manually with `stop_ortools!()`.

Use `solver="ORTools"` for minimum-cost circulation
or `max_flow_solver="ORTools"` for maximum-flow computations.

OR-Tools supports the Ford-Fulkerson algorithm (`"FF"`) and the integer-valued parametric search
methods (`"BS"`, `"IBS"`, `"INT"`, and `"INT-SG"`).

Newton-type search (`"NEWT"`) and the LP-based approach (`"LP"`) require HiGHS or Gurobi.


### Optional Gurobi Support

Users with a Gurobi license can install and load Gurobi:

```julia-repl
using Pkg
Pkg.add("Gurobi")

using DiFloT, Gurobi
```

Gurobi can then be selected via `solver="Gurobi"`.

## Quick Start

```julia
using DiFloT

# 1. Define a network via arrays
starts = [1, 1, 2, 3]
ends   = [2, 3, 4, 4]
caps   = [10, 5, 10, 10]
tts    = [1, 2, 1, 1]

net = NetworkInstance(starts, ends, caps, tts;
                      source=1, sink=4)

# 2. Solve the MDFP (T = 10)
mdf = solve_max_dynamic_flow(net, 1, 4, 10;
                             method="FF")

println("Maximum dynamic flow: ", mdf.value)

# 3. Solve the QFP (F = 20)
qf = solve_quickest_flow(net, 1, 4, 20;
                         method="IBS")

println("Minimum time-horizon: ", qf.time_T)
println("MCCP oracle calls: ", qf.mccp_calls)

@assert mdf.status == STATUS_OPTIMAL
@assert qf.status == STATUS_OPTIMAL
@assert mdf.value == 130.0
@assert qf.time_T == 3.0
```


## Solution Methods

### Maximum Dynamic Flow Problem

| Method | Description |
|--------|-------------|
| `"FF"` | Ford–Fulkerson algorithm |
| `"TE"` | Time-expanded network approach |

### Quickest Flow Problem

| Method | Description |
|--------|-------------|
| `"BS"` | Standard binary search |
| `"IBS"` | Improved binary search |
| `"INT"` | Interpolation search |
| `"INT-SG"` | Safeguarded interpolation search |
| `"NEWT"` | Newton-type search |
| `"LP"` | Lin–Jaillet LP-based approach |

The default methods are `"FF"` for the MDFP
and `"INT-SG"` for the QFP.

## Solver Configuration

The following example illustrates how to select the
optimization method and solver backend explicitly.

```julia
using DiFloT

net = NetworkInstance([1, 2], [2, 3],
                      [10, 10], [1, 1];
                      source=1, sink=3)

# Ford–Fulkerson with HiGHS dual simplex
mdf = solve_max_dynamic_flow(net, 1, 3, 10;
    method="FF",
    solver="HiGHS",
    solver_method="dual",
    with_paths=false)

# Safeguarded interpolation with HiGHS primal simplex
qf = solve_quickest_flow(net, 1, 3, 42;
    method="INT-SG",
    solver="HiGHS",
    solver_method="primal",
    with_paths=false)

@assert mdf.value == 90.0
@assert qf.time_T == 6.0
```

The option `with_paths=false` excludes path decomposition
from the computation.

For the QFP, the result additionally contains the number
of MCCP oracle calls (`mccp_calls`).
For safeguarded interpolation search, the number of
safeguard activations (`safeguard_activations`) is available.


## File Formats

`DiFloT.jl` uses a DIMACS-style format, separated by file extension for the two problem classes:
- `.min`: Standard DIMACS format.
- `.mdfp`: Extension to Maximum Dynamic Flow instances. Includes the time horizon `T` in the `p` line.
- `.qfp`: Extension to Quickest Flow instances. Includes the flow value `F` in the `p` line.

The `p` line defines the problem type, dimensions and problem parameter (`T`, the time-horizon for the MDFP, and `F`, the flow value for the QFP). The `n` lines define the source (`1`) and sink (`-1`). The `a` lines represent `(tail, head, lower_bound, capacity, transit_time)`.

### MDFP Instance
```text
c This is an MDFP instance
p MDFP 3 3 7
n 1 1
n 3 -1
a 1 2 0 10 1
a 2 3 0 10 1
a 1 3 0 5 3
```

### QFP Instance
```text
c This is a QFP instance
p QFP 3 3 42
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
using DiFloT

# Overrides the time horizon in the .mdfp file with T=50
res = solve_max_dynamic_flow_from_file(
    "data/network.mdfp"; s=1, t=3, T=50)

# Solves a .min file by providing T manually
res_min = solve_max_dynamic_flow_from_file(
    "data/network.min"; s=1, t=3, T=20)
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

## Running the Examples

The repository includes example scripts and network
instances in the `examples/` and `data/` directories.

From the repository root:

```bash
# Quickstart tutorial
julia --project=. examples/quickstart.jl

# Generate MDFP and QFP instances
julia --project=. examples/generate_instances.jl

# Illustrative benchmarks
julia --project=. examples/benchmark_mdfp.jl
julia --project=. examples/benchmark_qfp.jl
```


## Citing

If you use `DiFloT.jl` in your research, please cite the software repository.

A companion scientific publication is in preparation.

## Resources

- Documentation: [https://FelixBroesamle.github.io/DiFloT.jl/](https://FelixBroesamle.github.io/DiFloT.jl/).
- Package: [https://github.com/FelixBroesamle/DiFloT.jl](https://github.com/FelixBroesamle/DiFloT.jl).


## License

Distributed under the MIT License. See `LICENSE` for more information.

