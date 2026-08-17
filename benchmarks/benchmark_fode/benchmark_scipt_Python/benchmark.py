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

# Right-hand side of a 3D multi-order fractional system.
# The nonlinear term in the first equation is constructed so that
# a known analytic trajectory is available for verification.
def f(t: float, y: np.array) -> np.array:
    return np.array([
        1/math.sqrt(math.pi)*(((y[1]-0.5)*(y[2]-0.3))**(1/6) + math.sqrt(t)),
        math.gamma(2.2)*(y[0]-1),
        math.gamma(2.8)/math.gamma(2.2)*(y[1]-0.5),
    ])

# Analytic reference solution used to compute benchmark error.
# The second argument is unused and kept for callback compatibility.
def analytic(t: float, y: np.array) -> np.array:
    return np.array([t+1, t**1.2+0.5, t**1.8+0.3])

# Fractional derivative orders for each system component.
alpha = (0.5, 0.2, 0.6)

# Initial condition y(0).
y0 = np.array([1.0, 0.5, 0.3])

# Optional containers for trajectory storage (currently unused).
ts = []
ys = []

# Step-size sweep: dt = 2^{-i}, i = 3,...,7.
dts = [2.0**(-i) for i in range(3, 9)]

# Benchmark table with one row per dt:
#   - time: execution time in seconds
#   - error: norm of numerical minus analytic trajectory
df = pd.DataFrame({'time': [],
                   'error': []})

# Run benchmark for each step size.
for dt in dts:
    # Configure a fixed-step PECE Caputo solver.
    stepper = caputo.PECE(
        derivative_order=alpha,
        control=make_fixed_controller(dt, tstart=0.0, tfinal=5.0),
        source=partial(f),
        y0=(y0,),
        corrector_iterations=1,
    )

    # Measure wall-clock runtime of one complete solve.
    start_time = time.time()
    # Execute fractional evolution on [0, 5].
    solution = fracevolve(stepper, dtinit=dt)
    end_time = time.time()

    exec_time = end_time - start_time

    # Evaluate analytic solution on the solver time grid and compute error.
    ana = analytic_solution = analytic(solution.t, None)
    error = np.linalg.norm(solution.y - ana)

    # Append benchmark record for this dt.
    new_row = pd.DataFrame({'time': [exec_time],
                            'error': [error]})
    df = pd.concat([df, new_row], ignore_index=True)

# Export benchmark results to CSV for later plotting/comparison.
df.to_csv('/Users/quqingyu/SciFracX/paper/benchmarks/data/PYCAPUTO_PECE.csv')