#!/usr/bin/env python3
"""Run one pycaputo PECE solve for the cross-language FODE RSS benchmark.

The benchmark is the 10-state system D_t^0.8 y = A y on [0, 20], where A
has -1.5 on the diagonal and 0.25 on its first off-diagonals. The initial
state is [1.0, 0.9, ..., 0.1].

Usage:
    python3 memory_fode/fode_memory.py N
    python3 memory_fode/fode_memory.py --baseline

The final stdout line is parsed by ``memory_figure.jl``.  This worker times
only the numerical solve.  Its peak RSS is measured externally so that the
reported value includes Python, NumPy, pycaputo, and their loaded libraries.
"""

from __future__ import annotations

import gc
import os
import sys
import time

# Keep implicit numerical-library threading consistent with the other workers.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from pycaputo.controller import FixedController
from pycaputo.derivatives import CaputoDerivative
from pycaputo.events import StepCompleted
from pycaputo.fode import caputo
from pycaputo.stepping import evolve


ALPHA = 0.8
T0 = 0.0
TFINAL = 20.0
SYSTEM_SIZE = 10
Y0 = np.linspace(1.0, 0.1, SYSTEM_SIZE, dtype=float)
SYSTEM_MATRIX = (
    -1.5 * np.eye(SYSTEM_SIZE, dtype=float)
    + 0.25 * np.eye(SYSTEM_SIZE, k=1, dtype=float)
    + 0.25 * np.eye(SYSTEM_SIZE, k=-1, dtype=float)
)


def source(_t: float, y: np.ndarray) -> np.ndarray:
    """Right-hand side of the coupled linear fractional system."""

    return SYSTEM_MATRIX @ y


def solve_pece(
    number_of_steps: int, *, final_time: float = TFINAL
) -> tuple[np.ndarray, np.ndarray, float]:
    """Solve on the fixed interval using exactly the requested nominal h."""

    h = (final_time - T0) / number_of_steps
    stepper = caputo.PECE(
        ds=tuple(CaputoDerivative(ALPHA) for _ in range(SYSTEM_SIZE)),
        # Stop by the exact step count. When ``tfinal`` is supplied,
        # pycaputo's FixedController adds a small epsilon to every later step,
        # which accumulates into a measurable final-time overshoot for large N.
        control=FixedController(
            tstart=T0,
            tfinal=None,
            nsteps=number_of_steps,
            dt=h,
        ),
        source=source,
        y0=(Y0.copy(),),
        corrector_iterations=1,
    )

    times: list[float] = []
    values: list[np.ndarray] = []
    for event in evolve(stepper, dtinit=h):
        if isinstance(event, StepCompleted):
            times.append(float(event.t))
            values.append(np.array(event.y, copy=True))

    if not values:
        raise RuntimeError("pycaputo PECE returned no completed steps")

    times_array = np.asarray(times)
    values_array = np.stack(values)
    if times_array.size != number_of_steps + 1:
        raise RuntimeError("pycaputo PECE returned an unexpected number of time points")
    if values_array.shape != (number_of_steps + 1, SYSTEM_SIZE):
        raise RuntimeError("pycaputo PECE returned an unexpected state shape")
    time_tolerance = (
        8
        * number_of_steps
        * np.finfo(float).eps
        * max(1.0, abs(T0), abs(final_time))
    )
    if not np.isclose(times_array[-1], final_time, rtol=0.0, atol=time_tolerance):
        raise RuntimeError("pycaputo PECE did not reach the requested final time")

    return times_array, values_array, h


def main(argv: list[str]) -> None:
    if len(argv) != 2:
        raise SystemExit("usage: fode_memory.py N | --baseline")

    if argv[1] == "--baseline":
        # Imports above initialize the Python numerical stack.  Deliberately
        # skip the solver so the parent can measure its process-RSS baseline.
        print("FODE_MEMORY_RESULT,0,0,0", flush=True)
        return

    try:
        number_of_steps = int(argv[1])
    except ValueError as exc:
        raise SystemExit(f"N must be a positive integer; received {argv[1]!r}") from exc
    if number_of_steps <= 0:
        raise SystemExit("N must be positive")

    # Exercise imports and solver setup before timing.  The warm-up remains in
    # the measured process, which is appropriate for complete-process RSS.
    h = (TFINAL - T0) / number_of_steps
    warmup_steps = min(number_of_steps, 64)
    warmup_final_time = T0 + warmup_steps * h
    warmup_times, warmup_values, _ = solve_pece(
        warmup_steps, final_time=warmup_final_time
    )
    del warmup_times, warmup_values
    gc.collect()

    start = time.perf_counter()
    times, values, h = solve_pece(number_of_steps)
    runtime_seconds = time.perf_counter() - start

    if not np.all(np.isfinite(values[-1])):
        raise RuntimeError("pycaputo PECE returned a non-finite solution")

    # Keep the complete trajectory alive until after the timing/RSS sample.
    _trajectory_size = times.nbytes + values.nbytes
    print(
        f"FODE_MEMORY_RESULT,{number_of_steps},{h:.17g},{runtime_seconds:.17g}",
        flush=True,
    )
    del _trajectory_size


if __name__ == "__main__":
    main(sys.argv)
