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
from julia import FractionalDiffEq

def f(t: float, y: np.array) -> np.array:
    return np.array([
            -10*y[0]
    ])
alpha = 0.8
def analytic(t: float, y: np.array) -> np.array:
    return np.array([FractionalDiffEq.mittleff(alpha, -10*t**alpha)])


y0 = np.array([1.0])

ts = []
ys = []

dts = [2.0**(-i) for i in range(3, 8)]
df = pd.DataFrame({'time': [],
                   'error': []})
for dt in dts:
    stepper = caputo.PECE(
        derivative_order=alpha,
        control=make_fixed_controller(dt, tstart=0.0, tfinal=5.0),
        source=partial(f),
        y0=(y0,),
        corrector_iterations=1,
    )
    start_time = time.time()
    # Code to benchmark
    solution = fracevolve(stepper, dtinit=dt)
    end_time = time.time()

    exec_time = end_time - start_time

    ana = analytic(solution.t, None)
    error = np.linalg.norm(solution.y - ana)
    new_row = pd.DataFrame({'time': [exec_time],
                            'error': [error]})
    df = pd.concat([df, new_row], ignore_index=True)

df.to_csv('/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Singleterm_PYCAPUTO_PECE.csv')