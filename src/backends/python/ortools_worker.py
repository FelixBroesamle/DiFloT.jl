"""Private binary IPC worker; never imported into Julia's native-library process."""
import sys
import struct
import traceback
import os

# Reserve stdout for the binary protocol; suppress import diagnostics.
inp, out = sys.stdin.buffer, sys.stdout.buffer
sys.stdout = open(os.devnull, "w")

import numpy as np
import ortools
from ortools.graph.python import max_flow, min_cost_flow

# Little-endian int64 request and response headers.
HEADER = struct.Struct('<5q')
REPLY = struct.Struct('<3q')


def read_exact(size):
    data = bytearray(size)
    view = memoryview(data)
    pos = 0
    while pos < size:
        count = inp.readinto(view[pos:])
        
        if not count:
            raise EOFError('Truncated request')
        pos += count

    return data


def array(count):
    return np.frombuffer(read_exact(8 * count), dtype='<i8')


def solve(op, n, m, s, t, tails, heads, caps, costs=None):
    if not (0 <= s < n and 0 <= t < n and s != t and n <= 2**31-1):
        raise ValueError('Invalid terminal or node count')
    
    if (np.any(tails < 0) or np.any(heads < 0)
            or np.any(tails >= n) or np.any(heads >= n)
            or np.any(caps < 0)):
        raise ValueError('Invalid arc data')
    
    # OR-Tools uses int32 node indices and int64 capacities/costs.
    tails, heads = tails.astype(np.int32), heads.astype(np.int32)

    # Operation 1: maximum flow.
    if op == 1:
        solver = max_flow.SimpleMaxFlow()
        solver.add_arcs_with_capacity(tails, heads, caps)
        # Zero-capacity terminal arcs ensure even isolated terminals exist.
        solver.add_arc_with_capacity(s, s, 0)
        solver.add_arc_with_capacity(t, t, 0)
        status = solver.solve(s, t)
    
        if status != solver.OPTIMAL:
            raise RuntimeError(f'SimpleMaxFlow status: {status}')
    
        return int(solver.optimal_flow()), np.empty(0, dtype='<i8')

    # Operation 2: minimum-cost circulation.
    if op == 2:
        solver = min_cost_flow.SimpleMinCostFlow()
        indices = solver.add_arcs_with_capacity_and_unit_cost(tails, heads, caps, costs)
    
        # All supplies are zero by default (circulation). Explicitly create last node.
        solver.set_node_supply(n - 1, 0)
        status = solver.solve()
    
        if status != solver.OPTIMAL:
            raise RuntimeError(f'SimpleMinCostFlow status: {status}')
    
        return -int(solver.optimal_cost()), np.asarray(solver.flows(indices), dtype='<i8')
    
    raise ValueError(f'Unknown operation {op}')


if __name__ == '__main__':
    out.write(f'DIFLOT-ORTOOLS/1 {ortools.__version__}\n'.encode('utf-8'))
    out.flush()

    # Serve requests until Julia closes the stream or sends operation 0.
    while True:
        try:
            header = inp.read(HEADER.size)
    
            if not header:
                break
            if len(header) != HEADER.size:
                raise EOFError('Truncated header')
    
            op, n, m, s, t = HEADER.unpack(header)
            if op == 0:
                break
            if op not in (1, 2) or m < 0:
                raise ValueError('Invalid request header')
    
            tails, heads, caps = array(m), array(m), array(m)
            costs = array(m) if op == 2 else None
            value, flows = solve(op, n, m, s, t, tails, heads, caps, costs)
    
            out.write(REPLY.pack(0, value, len(flows)))
            out.write(flows.tobytes())
            out.flush()
    
            del tails, heads, caps, costs, flows
    
        except Exception as exc:
            traceback.print_exc(file=sys.stderr)
            msg = str(exc).encode('utf-8', errors='replace')
    
            out.write(REPLY.pack(1, 0, len(msg)))
            out.write(msg)
            out.flush()
    
            break  # Do not try to recover from an unknown protocol position.
