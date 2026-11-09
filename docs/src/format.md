# The DIMACS .min, .mdfp, and .qfp Data Formats

`DiFloT.jl` uses a structured text standard based on the classic DIMACS minimum-cost flow format (`.min`) to natively capture discrete-time dynamic flow problems. The format is separated by file extension to distinguish between the two problem classes.

## File Structure Overview

A standard file consists of three types of lines:
* **Comment Line (`c`)**: Optional comments.
* **Problem Line (`p`)**: Defines the problem type and dimensions.
* **Node Line (`n`)**: Defines source (`1`) and sink (`-1`).
* **Arc Line (`a`)**: Defines the arcs `(tail, head, lower_bound, capacity, transit_time)`.

## 1. The Base Format (`.min`)
The `.min` format defines the static topology of the network.
```text
c This is a standard DIMACS .min instance
p min 3 3
n 1 1
n 3 -1
a 1 2 0 10 1
a 2 3 0 10 1
a 1 3 0 5 3
```

## 2. The Maximum Dynamic Flow Format (.mdfp)

Extends .min by appending the time horizon `T` to the problem line.

```text
c This is an MDFP instance
p MDFP 3 3 10
n 1 1
n 3 -1
a 1 2 0 10 1
a 2 3 0 10 1
a 1 3 0 5 3
```

The `p` line format is: `p MDFP <num_nodes> <num_arcs> <time_horizon_T>

### 3. The Quickest Flow Format (.qfp)

Extends `.min` by appending the target flow value `F` to the problem line.

```text
c This is a QFP instance
p QFP 3 3 20
n 1 1
n 3 -1
a 1 2 0 10 1
a 2 3 0 10 1
a 1 3 0 5 3
```

The `p` line format is: `p QFP <num_nodes> <num_arcs> <target_flow_F>

##






