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
from optimizationParameters import numberOfOptimizations

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
    global_sol = None
    for i in range(numberOfAttempts):
        newGuess = getInitialGuess()
        try:
            local_sol = optimize(timeseries=timeseries, x0=newGuess)
        except Exception:
            print("An error occured during the optimization with the following parameters")
            printParameters(newGuess)
            print("Try new initial guess")
            continue
        print("local minimum with "+str(local_sol.fun))
        print(str(i)+" th call of optimize")
        printParameters(newGuess)
        print("----------------------------------------")
        if global_sol is None or local_sol.fun < global_sol.fun:
            global_sol = local_sol
            initialGuess = newGuess
    return (global_sol, initialGuess)

def printParameters(parameters):
    """
    Prints parameters in internal and human readable format
    :param parameters:
    :return: -
    """
    print("Parameters: "+str(parameters))
    k, theta, d = convertToHumanReadableParameters(parameters=parameters)
    print("In human readable format")
    print("k_initial: "+str(k))
    print("theta_initial: "+str(theta))
    print("d_initial: "+str(d))



sol, initialGuess = optimizeSeveralTimes(timeseries, numberOfOptimizations)
print("--------------------------------------------------------------")
print("Result of computation:")
print(sol)
print("With initial guess")
printParameters(initialGuess)


(time, initialValue) = timeseries[0]
(stoppingTime, p) = timeseries[-1]
plotODEsolution(parameters=sol.x, initialValue=initialValue, stoppingTime=stoppingTime, numberOfSteps=numberOfSteps,
                method=methodForODE)

