#!/usr/bin/env python3
"""Run one pycaputo PECE dimension/interval scalability benchmark.

Usage:
    python3 scalability/fode_scalability.py D T N
    python3 scalability/fode_scalability.py --baseline

The benchmark is D_t^0.8 y = A y with the same matrix, initial state, time
grid, and one PECE corrector iteration used by the Julia and MATLAB workers.
"""

from __future__ import annotations

import gc
import os
import sys
import time

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


def system_matrix(dimension: int) -> np.ndarray:
    return (
        -1.5 * np.eye(dimension, dtype=float)
        + 0.25 * np.eye(dimension, k=1, dtype=float)
        + 0.25 * np.eye(dimension, k=-1, dtype=float)
    )


def initial_state(dimension: int) -> np.ndarray:
    return np.linspace(1.0, 0.1, dimension, dtype=float)


def solve_pece(
    dimension: int, final_time: float, number_of_steps: int
) -> tuple[np.ndarray, np.ndarray, float]:
    h = (final_time - T0) / number_of_steps
    matrix = system_matrix(dimension)

    def source(_t: float, y: np.ndarray) -> np.ndarray:
        return matrix @ y

    stepper = caputo.PECE(
        ds=tuple(CaputoDerivative(ALPHA) for _ in range(dimension)),
        # Stop by N rather than tfinal. Supplying tfinal makes pycaputo add a
        # small epsilon to later steps, so the resulting mesh is not exactly h.
        control=FixedController(
            tstart=T0,
            tfinal=None,
            nsteps=number_of_steps,
            dt=h,
        ),
        source=source,
        y0=(initial_state(dimension),),
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
        raise RuntimeError("pycaputo PECE returned an unexpected number of points")
    if values_array.shape != (number_of_steps + 1, dimension):
        raise RuntimeError("pycaputo PECE returned an unexpected state shape")
    tolerance = (
        8
        * number_of_steps
        * np.finfo(float).eps
        * max(1.0, abs(T0), abs(final_time))
    )
    if not np.isclose(times_array[-1], final_time, rtol=0.0, atol=tolerance):
        raise RuntimeError("pycaputo PECE did not reach the requested final time")

    return times_array, values_array, h


def parse_positive_int(value: str, name: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise SystemExit(f"{name} must be a positive integer") from exc
    if parsed <= 0:
        raise SystemExit(f"{name} must be positive")
    return parsed


def parse_positive_float(value: str, name: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise SystemExit(f"{name} must be a positive number") from exc
    if not np.isfinite(parsed) or parsed <= T0:
        raise SystemExit(f"{name} must be finite and positive")
    return parsed


def main(argv: list[str]) -> None:
    if argv[1:] == ["--baseline"]:
        print("SCALABILITY_RESULT,0,0,0,0,0", flush=True)
        return
    if len(argv) != 4:
        raise SystemExit("usage: fode_scalability.py D T N | --baseline")

    dimension = parse_positive_int(argv[1], "D")
    final_time = parse_positive_float(argv[2], "T")
    number_of_steps = parse_positive_int(argv[3], "N")
    h = (final_time - T0) / number_of_steps

    warmup_steps = min(number_of_steps, 64)
    warmup_final_time = T0 + warmup_steps * h
    warmup_times, warmup_values, _ = solve_pece(
        dimension, warmup_final_time, warmup_steps
    )
    del warmup_times, warmup_values
    gc.collect()

    start = time.perf_counter()
    times, values, returned_h = solve_pece(dimension, final_time, number_of_steps)
    runtime_seconds = time.perf_counter() - start
    if returned_h != h:
        raise RuntimeError("internal step-size mismatch")
    if not np.all(np.isfinite(values[-1])):
        raise RuntimeError("pycaputo PECE returned a non-finite solution")

    # Retain the returned trajectory until after the RSS sample/result output.
    _trajectory_size = times.nbytes + values.nbytes
    print(
        "SCALABILITY_RESULT,"
        f"{dimension},{final_time:.17g},{number_of_steps},"
        f"{h:.17g},{runtime_seconds:.17g}",
        flush=True,
    )
    del _trajectory_size


if __name__ == "__main__":
    main(sys.argv)
