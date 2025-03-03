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

def f(t: float, y: np.array) -> np.array:
    return np.array([
        1/math.sqrt(math.pi)*(((y[1]-0.5)*(y[2]-0.3))**(1/6) + math.sqrt(t)),
        math.gamma(2.2)*(y[0]-1),
        math.gamma(2.8)/math.gamma(2.2)*(y[1]-0.5),
    ])

def analytic(t: float, y: np.array) -> np.array:
    return np.array([t+1, t**1.2+0.5, t**1.8+0.3])

alpha = (0.5, 0.2, 0.6)

y0 = np.array([1.0, 0.5, 0.3])

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

    ana = analytic_solution = analytic(solution.t, None)
    error = np.linalg.norm(solution.y - ana) #TODO: see how the norm is working?
    new_row = pd.DataFrame({'time': [exec_time],
                            'error': [error]})
    df = pd.concat([df, new_row], ignore_index=True)

df.to_csv('/Users/quqingyu/SciFracX/paper/benchmarks/data/PYCAPUTO_PECE.csv')