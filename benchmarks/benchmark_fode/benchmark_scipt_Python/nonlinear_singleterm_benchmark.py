from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from pycaputo import fracevolve, fracplot
from pycaputo.controller import make_fixed_controller
from pycaputo.fode import caputo
from pycaputo.events import StepCompleted
from pycaputo.stepping import evolve
import time
import math
import pandas as pd

# Nonlinear right-hand side for a single-term Caputo FODE.
# The forcing is manufactured so that an analytic solution is available.
def f(t: float, y: np.array) -> np.array:
    return np.array([
            (40320 / math.gamma(9 - 0.5) * t ** (8 - 0.5) - 3 * math.gamma(5 + 0.5 / 2)/ math.gamma(5 - 0.5 / 2) * t ** (4 - 0.5 / 2) + 9/4 * math.gamma(0.5 + 1) +(3 / 2 * t ** (0.5 / 2) - t ** 4) ** 3 - y[0] ** (3 / 2))
    ])

# Fractional derivative order.
alpha = 0.5

# Closed-form reference solution used for error evaluation.
# Signature includes y to match generic callback style, though y is unused.
def analytic(t: float, y: np.array) -> np.array:
    return np.array([t**8 - 3 * t ** (4 + alpha / 2) + 9 / 4 * t**alpha])


# Initial condition y(0).
y0 = np.array([0.0])

# Auxiliary containers (not used later, kept for compatibility/extension).
ts = []
ys = []

# Benchmark step sizes: dt = 2^{-i}, i = 3,...,7.
dts = [2.0**(-i) for i in range(3, 8)]

# Output table with one row per dt: execution time and global error norm.
df = pd.DataFrame({'time': [],
                   'error': []})

# Sweep all step sizes and collect runtime/accuracy statistics.
for dt in dts:
    # Configure fixed-step PECE solver for Caputo derivative.
    stepper = caputo.PECE(
        derivative_order=alpha,
        control=make_fixed_controller(dt, tstart=0.0, tfinal=1.0),
        source=partial(f),
        y0=(y0,),
        corrector_iterations=1,
    )

    # Measure wall-clock execution time for a single full solve.
    start_time = time.time()
    # Run the fractional evolution.
    solution = fracevolve(stepper, dtinit=dt)
    end_time = time.time()

    exec_time = end_time - start_time

    # Compute error against analytic solution on the solver's time grid.
    ana = analytic(solution.t, None)
    error = np.linalg.norm(solution.y - ana)

    # Append one benchmark record for the current dt.
    new_row = pd.DataFrame({'time': [exec_time],
                            'error': [error]})
    df = pd.concat([df, new_row], ignore_index=True)

# Export benchmark results to CSV.
# Columns: time (seconds), error (L2 norm over trajectory samples).
df.to_csv('/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_PYCAPUTO_PECE.csv')