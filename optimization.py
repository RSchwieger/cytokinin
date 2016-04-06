__author__ = 'Robert Schwieger'

from generalODESystem import computeError
from generalODESystem import plotODEsolution
from scipy.optimize import root
from scipy.optimize import minimize


#------------
from optimizationParameters import getBounds
from optimizationParameters import getInitialGuess as getInitialGuess
#------------

from optimizationParameters import method
from optimizationParameters import numberOfSteps
from optimizationParameters import methodForODE
from optimizationParameters import timeseries

def optimize(timeseries):
    (time, p) = timeseries[-1]
    dt = time/numberOfSteps
    def func(parameters):
        return computeError(timeseries=timeseries, parameters=parameters, dt=dt, method=methodForODE)
    bounds = getBounds()
    sol = minimize(fun=func, x0=getInitialGuess(), method=method, bounds=bounds)
    # sol = root(fun=func, x0=getInitialGuess(), jac=False, method='lm')
    return sol


sol = optimize(timeseries)
print(sol)


