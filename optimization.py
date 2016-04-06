__author__ = 'Robert Schwieger'

from generalODESystem import computeError
from generalODESystem import plotODEsolution
from scipy.optimize import root
from scipy.optimize import minimize


#------------
from optimizationParameters import getBounds
from optimizationParameters import getInitialGuess
from optimizationParameters import convertToHumanReadableParameters
#------------

from optimizationParameters import method
from optimizationParameters import numberOfSteps
from optimizationParameters import methodForODE
from optimizationParameters import timeseries
from optimizationParameters import x0
from optimizationParameters import ftol

def optimize(timeseries, x0):
    """
    Finds the best curve fit in the set of ODE-solutions to a specified time series data.
    :param timeseries: time series data
    :param x0: initial parameter guess
    :return:OptimizeResult
    """
    (time, p) = timeseries[-1]
    dt = time/numberOfSteps
    def func(parameters):
        return computeError(timeseries=timeseries, parameters=parameters, dt=dt, method=methodForODE)
    bounds = getBounds()
    sol = minimize(fun=func, x0=x0, method=method, bounds=bounds, options={'ftol':ftol})
    return sol

def optimizeSeveralTimes(timeseries, numberOfAttempts):
    """
    Wrapper to run optimize several times to increase the chance to find a global minimum
    :param timeseries: time series data
    :param numberOfAttempts: number of optimize-calls
    :return: tuple: (best solution, corresponding initialGuess)
    """
    initialGuess = getInitialGuess()
    global_sol = optimize(timeseries=timeseries, x0=initialGuess)
    for i in range(numberOfAttempts-1):
        newGuess = getInitialGuess()
        local_sol = optimize(timeseries=timeseries, x0=newGuess)
        print("local minimum with "+str(local_sol.fun))
        if local_sol.fun < global_sol.fun:
            global_sol = local_sol
            initialGuess = newGuess
    return (global_sol, initialGuess)


sol, initialGuess = optimizeSeveralTimes(timeseries, 5)
print(sol)
print("With initial guess: "+str(initialGuess))
k, theta, d = convertToHumanReadableParameters(x0)
print("In human readable format")
print("k_initial: "+str(k))
print("theta_initial: "+str(theta))
print("d_initial: "+str(d))

(time, initialValue) = timeseries[0]
(stoppingTime, p) = timeseries[-1]
plotODEsolution(parameters=sol.x, initialValue=initialValue, stoppingTime=stoppingTime, numberOfSteps=numberOfSteps,
                method=methodForODE)

