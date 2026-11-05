# DiFloT.jl

**DiFloT.jl** (**Di**screte **Flo**ws over **T**ime) is a high-performance Julia package for solving the discrete-time Maximum Dynamic Flow Problem (MDFP) and Quickest Flow Problem (QFP).

## Key Features

It implements the following solution methods:
1. For the MDFP:
    - The Ford-Fulkerson Algorithm (`"FF"`).
    - Solving the Problem on the time-expanded network (`"TE"`).
2. For the QFP:
    - A standard binary search (`"BS"`).
    - Improved binary search (`"IBS"`).
    - Interpolation search (`"INT"`) and a safeguarded Interpolation search (`"INT-SG"`).
    - Newton search (`"NEWT"`).
    - Lin-Jaillet LP formulation (`"LP"`).

## Installation

### Standard installation
```julia
using Pkg
Pkg.add("DiFloT")
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

### Running the Example Workflows

See the `examples/` directory for Jupyter notebooks and scripts.


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
mdf = solve_max_dynamic_flow(net, 1, 4, 10; method="FF", solver="HiGHS")
println("Max Flow: ", mdf.value)

# 3. Solve Quickest Flow Problem (QFP)
qf = solve_quickest_flow(net, 1, 4, 20; method="IBS", solver="HiGHS")
println("Quickest Time: ", qf.time_T)
```

## File Format

`DiFloT.jl` uses a DIMACS-style format, separated by file extension for the two problem classes:
- `.mdfp` for Maximum Dynamic Flow instances.
- `.qfp` for Quickest Flow instances.

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

### From File

```julia
net = parse_min("network.min")
res = solve_max_dynamic_flow(net, s=1, t=3, T=10)
```

### In-Memory Generation

You can define networks directly using arrays:

```julia
start_nodes = [1, 1, 2, 2, 3]
end_nodes   = [2, 3, 3, 4, 4]
capacities  = [15, 8, 20, 4, 10]
transit_times = [4, 4, 2, 2, 6]

net = NetworkInstance(4, start_nodes, end_nodes, capacities, transit_times)
res = solve_quickest_flow(net, 1, 4, 50.0, method="IBS")
```

## Citing

If you use `DiFloT.jl` in your research, please cite the following paper for the methodologies:

```text
To appear.
```

And the following paper for citing the Software:
```text
To appear.
```

## Resources

- Documentation: (To appear). [https://FelixBroesamle.github.io/DiFloT.jl/](https://FelixBroesamle.github.io/DiFloT.jl/).
- Package: [https://github.com/FelixBroesamle/DiFloT.jl](https://github.com/FelixBroesamle/DiFloT.jl).

## License

Distributed under the MIT License. See `LICENSE` for more information.








