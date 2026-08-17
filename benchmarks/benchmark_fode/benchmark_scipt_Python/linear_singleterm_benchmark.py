import time
import numpy as np
import pandas as pd

from pycaputo.controller import make_fixed_controller
from pycaputo.fode import caputo
from pycaputo.derivatives import CaputoDerivative as D
from pycaputo.events import StepCompleted
from pycaputo.stepping import evolve

from julia.api import Julia

# ------------------------------------------------------------
# Initialize Julia
# ------------------------------------------------------------
jl = Julia(compiled_modules=False)

from julia import FractionalDiffEq


# ------------------------------------------------------------
# Problem:
#
#   D_t^alpha y(t) = -10 y(t)
#   y(0) = 1
#
# Exact solution:
#
#   y(t) = E_alpha(-10 t^alpha)
# ------------------------------------------------------------

alpha = 0.8
lambda_ = -10.0


def f(t: float, y: np.ndarray) -> np.ndarray:
    return np.array([
        lambda_ * y[0]
    ])


y0 = np.array([1.0])


# ------------------------------------------------------------
# Analytic solution
# ------------------------------------------------------------
def analytic(t):
    """
    Evaluate the analytical solution

        y(t) = E_alpha(-10 t^alpha)

    using FractionalDiffEq.jl.
    """

    t = np.asarray(t, dtype=float)

    y_exact = np.array([
        float(
            FractionalDiffEq.mittleff(
                alpha,
                lambda_ * float(ti) ** alpha
            )
        )
        for ti in t
    ])

    return y_exact


# ------------------------------------------------------------
# Step sizes
# ------------------------------------------------------------
dts = [2.0 ** (-i) for i in range(3, 9)]


# ------------------------------------------------------------
# Benchmark results
# ------------------------------------------------------------
results = []


for dt in dts:

    # --------------------------------------------------------
    # Construct PECE solver
    # --------------------------------------------------------
    stepper = caputo.PECE(

        # IMPORTANT:
        # ds must be a tuple in the current pycaputo API.
        #
        # One equation -> one Caputo derivative
        ds=(D(alpha),),

        control=make_fixed_controller(
            dt,
            tstart=0.0,
            tfinal=5.0
        ),

        source=f,

        # Initial values are also stored as a tuple.
        # For alpha in (0, 1), only y(0) is needed.
        y0=(y0,),

        corrector_iterations=1,
    )

    # --------------------------------------------------------
    # Solve and benchmark
    # --------------------------------------------------------
    ts = []
    ys = []

    start_time = time.perf_counter()

    for event in evolve(stepper):

        # evolve can in principle generate several event types.
        # We only collect successfully completed time steps.
        if isinstance(event, StepCompleted):
            ts.append(event.t)
            ys.append(event.y.copy())

    end_time = time.perf_counter()

    exec_time = end_time - start_time

    # --------------------------------------------------------
    # Convert solution to NumPy arrays
    # --------------------------------------------------------
    ts = np.asarray(ts)

    # shape before squeeze:
    #     (N, 1)
    ys = np.asarray(ys)

    # scalar equation -> shape (N,)
    y_num = ys[:, 0]

    # --------------------------------------------------------
    # Analytical solution
    # --------------------------------------------------------
    y_exact = analytic(ts)

    # --------------------------------------------------------
    # Error
    # --------------------------------------------------------
    error_l2 = np.linalg.norm(y_num - y_exact)

    error_linf = np.max(np.abs(y_num - y_exact))

    results.append({
        "dt": dt,
        "N": len(ts),
        "time": exec_time,
        "error": error_linf,
    })


# ------------------------------------------------------------
# Save benchmark
# ------------------------------------------------------------
df = pd.DataFrame(results)

output_file = (
    "/Users/quqingyu/SciFracX/paper/benchmarks/data/"
    "Linear_Singleterm_PYCAPUTO_PECE.csv"
)

df.to_csv(
    output_file,
    index=False
)