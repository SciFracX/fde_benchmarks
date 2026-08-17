import time
import math
import numpy as np
import pandas as pd

from pycaputo.controller import make_fixed_controller
from pycaputo.derivatives import CaputoDerivative as D
from pycaputo.fode import caputo
from pycaputo.events import StepCompleted
from pycaputo.stepping import evolve


# ------------------------------------------------------------
# Fractional derivative order
# ------------------------------------------------------------
alpha = 0.5


# ------------------------------------------------------------
# Nonlinear single-term Caputo FODE
#
#   D_t^alpha y(t) = f(t, y)
#
# with manufactured exact solution
#
#   y(t) = t^8 - 3 t^(4 + alpha/2) + 9/4 t^alpha
#
# and
#
#   y(0) = 0
# ------------------------------------------------------------
def f(t: float, y: np.ndarray) -> np.ndarray:
    return np.array([
        (
            40320.0 / math.gamma(9.0 - alpha)
            * t ** (8.0 - alpha)

            - 3.0
            * math.gamma(5.0 + alpha / 2.0)
            / math.gamma(5.0 - alpha / 2.0)
            * t ** (4.0 - alpha / 2.0)

            + 9.0 / 4.0
            * math.gamma(alpha + 1.0)

            + (
                3.0 / 2.0 * t ** (alpha / 2.0)
                - t ** 4
            ) ** 3

            - y[0] ** (3.0 / 2.0)
        )
    ])


# ------------------------------------------------------------
# Exact solution
# ------------------------------------------------------------
def analytic(t):
    t = np.asarray(t, dtype=float)

    return (
        t ** 8
        - 3.0 * t ** (4.0 + alpha / 2.0)
        + 9.0 / 4.0 * t ** alpha
    )


# ------------------------------------------------------------
# Initial condition
# ------------------------------------------------------------
y0 = np.array([0.0])


# ------------------------------------------------------------
# Benchmark step sizes
#
# dt = 2^{-i}, i = 3,...,8
# ------------------------------------------------------------
dts = [2.0 ** (-i) for i in range(3, 9)]


# ------------------------------------------------------------
# Benchmark results
# ------------------------------------------------------------
results = []


for dt in dts:


    # --------------------------------------------------------
    # Construct PECE solver
    #
    # IMPORTANT:
    #
    # New pycaputo API:
    #
    #   ds=(D(alpha),)
    #
    # instead of:
    #
    #   derivative_order=alpha
    #
    # --------------------------------------------------------
    stepper = caputo.PECE(

        ds=(D(alpha),),

        control=make_fixed_controller(
            dt,
            tstart=0.0,
            tfinal=1.0,
        ),

        source=f,

        y0=(y0,),

        corrector_iterations=1,
    )


    # --------------------------------------------------------
    # Solve
    #
    # evolve(stepper) now returns events rather than a
    # solution object with .t and .y.
    # --------------------------------------------------------
    ts = []
    ys = []

    start_time = time.perf_counter()

    for event in evolve(stepper):

        if isinstance(event, StepCompleted):
            ts.append(event.t)
            ys.append(event.y.copy())

    end_time = time.perf_counter()

    exec_time = end_time - start_time


    # --------------------------------------------------------
    # Convert numerical solution to numpy arrays
    # --------------------------------------------------------
    ts = np.asarray(ts, dtype=float)
    ys = np.asarray(ys)

    # For this scalar FODE:
    #
    # ys.shape = (N, 1)
    #
    # Convert to:
    #
    # y_num.shape = (N,)
    #
    y_num = ys[:, 0]


    # --------------------------------------------------------
    # Exact solution
    # --------------------------------------------------------
    y_exact = analytic(ts)


    # --------------------------------------------------------
    # Error
    # --------------------------------------------------------

    # Same L2 norm definition as your original benchmark
    error_l2 = np.linalg.norm(
        y_num - y_exact
    )

    # Also calculate Linfinity for comparison
    error_linf = np.max(
        np.abs(y_num - y_exact)
    )

    results.append({
        "dt": dt,
        "N": len(ts),
        "time": exec_time,
        "error": error_linf,
    })


# ------------------------------------------------------------
# Save benchmark results
# ------------------------------------------------------------
df = pd.DataFrame(results)


output_file = (
    "/Users/quqingyu/SciFracX/paper/benchmarks/data/"
    "Singleterm_PYCAPUTO_PECE.csv"
)


df.to_csv(
    output_file,
    index=False,
)