__author__ = 'Robert Schwieger'

from generalODESystem import leastSquareError
from generalODESystem import plotODEsolution
from generalODESystem import computeSolutionOfODEAtTimepoints
from scipy.optimize import root
from scipy.optimize import minimize
from random import uniform
import matplotlib.pyplot as plt


#------------
from optimizationParameters import getBounds
from optimizationParameters import getInitialGuess
from optimizationParameters import convertToHumanReadableParameters
#------------

from optimizationParameters import method
from optimizationParameters import numberOfSteps
from optimizationParameters import methodForODE
from optimizationParameters import timeseries
from optimizationParameters import ftol
from optimizationParameters import numberOfOptimizations
from optimizationParameters import x0
from scipy.spatial import distance

from optimizationParameters import brute
from optimizationParameters import basinhopping

def objectiveFunction(timeseries, parameters, methodForODE, methodForOptimization):
    """
    Objective function of the minimization. Encapsulates a further minimization.
    :param timeseries: Given time series data for the fitting
    :param parameters: Values where to evaluate the objective function
    :param dt: precision
    :param method: Algorithm for minimization
    :return: objective value
    """
    # extract timepoints for comparison
    timePoints = [t for (t,p) in timeseries]
    # extract measurements
    dataPoints = [p for (t,p) in timeseries]
    # compute step size
    dt = timePoints[-1]/numberOfSteps
    # solve the normalized ODE-system with the given parameters
    ODETimeSeries, maxVector = computeSolutionOfODEAtTimepoints(initialValue=dataPoints[0], timePoints=timePoints, dt=dt,
                                                     parameters=parameters, method=methodForODE)
    # extract computed data for comparison with measurements
    computedPoints = [p for (t,p) in ODETimeSeries]

    # helper function
    def helperFunction(v):
        vPoints = [[v[i]*computedPoint[i] for i in range(len(v))] for computedPoint in computedPoints]
        return sum([distance.euclidean(dataPoints[i],vPoints[i])**2 for i in range(len(vPoints))])

    # compute boundaries based on the maximum values of the components of the solution
    bounds = [(0, 1/maxVector[i]) for i in range(len(dataPoints[0]))] #ToDo: zero Exception
    x0 = [uniform(a,b) for (a,b) in bounds] # Compute initial value for v randomly
    #x0 = [1 for (a,b) in bounds]
    # minimize the boundary problem

    v = [0 for x in x0] # initialize list of same length as x0 with arbitrary values
    for i in range(len(computedPoints[0])):
        a = [computedPoint[i] for computedPoint in computedPoints]
        b = [dataPoint[i] for dataPoint in dataPoints]
        def helperFunction1D(v):
            return sum([(v*a[i]-b[i])**2 for i in range(len(a))])
        sol = minimize(fun=helperFunction1D, x0=x0[i], method=methodForOptimization, bounds=[bounds[i]], options={'ftol':ftol})
        v[i] = sol.x[0]
    # print(helperFunction([1,1,1,1]))
    return helperFunction(v)

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
        # return objectiveFunction(timeseries=timeseries, parameters=parameters, methodForODE=methodForODE,
        #                  methodForOptimization=method)
        return leastSquareError(timeseries=timeseries, parameters=parameters, dt=dt, method=methodForODE)
    bounds = getBounds()
    sol = minimize(fun=func, x0=x0, method=method, bounds=bounds, options={'ftol':ftol})
    # sol = brute(func=func, ranges=tuple(bounds))
    return sol

def optimizeSeveralTimes(timeseries, numberOfAttempts):
    """
    Wrapper to run optimize several times to increase the chance to find a global minimum
    :param timeseries: time series data
    :param numberOfAttempts: number of optimize-calls
    :return: tuple: (best solution, corresponding initialGuess)
    """
    global_sol = None
    initialGuess = None
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
        print("Initial Guess: ")
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


parameters1 = [  6.18590237,   3.17409756,  21.73556717,   3.02734518,
         3.97093404,   0.45581197,   0.21500462,   0.02606325,
         0.05210057,   0.30661293,   0.19557138,   0.13829967,
         0.02742286,   0.24319961]
res = objectiveFunction(timeseries, parameters=parameters1,methodForODE=methodForODE, methodForOptimization=method)
print(res)

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