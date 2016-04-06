__author__ = 'Robert Schwieger'

from scipy.integrate import ode
from scipy.spatial import distance
import matplotlib.pyplot as plt
from optimizationParameters import convertToHumanReadableParameters, continuousHomolog

# Functions of the ODE-system


def f(t,x, parameters):
    """
    Function of the ODE-system written in the form required by the ODE-solver
    :param t: Not used
    :param x: list of input values
    :return: function evaluation saves as list
    """
    (k, theta, d) = convertToHumanReadableParameters(parameters=parameters)
    discreteTimestep = continuousHomolog(x, k, theta)
    return [d[i]*(discreteTimestep[i]-x[i]) for i in range(len(discreteTimestep))]


# Solving the ODE-system

def getSolutionAtTimeT(initialValue, stoppingTime, dt, parameters, method):
    """
    Solves the ODE-system and returns the state at the stopping time
    :param initialValue: initial Value of the ODE-system
    :param stoppingTime:
    :param dt: time step
    :return: state at stopping time
    """
    y = ode(f).set_integrator(method)
    y.set_initial_value(initialValue, 0.0).set_f_params(parameters) # set initial value at time = 0
    i = 0
    while y.successful() and y.t < stoppingTime:
        y.integrate(y.t+dt)
        if y.successful() is False:
            print("Something went wrong during r.integrate()")
            print("params: "+str(y.f_params))
            print("t: "+str(y.t))
            print("dt: "+str(dt))
            print("y.y: "+str(y.y))
            print("stoppingTime: "+str(stoppingTime))
            k, theta, d = convertToHumanReadableParameters(parameters)
            print(k)
            print(theta)
            print(d)
    return y

def getSolutionAtTimepoints(initialValue, timePoints, dt, parameters, method):
    """
    Computes a solution vector that contains only the values at the specified timepoints.
    :param initialValue: initial value of the ODE-System
    :param timePoints: monotonous list of time points
    :param dt: size of time steps
    :return: list of tuples. The first entry of a tuple is the time point, the second the corresponding value of the function
    """
    # compute a list of distances between the time steps
    relativeTimePoints = [timePoints[i+1]-timePoints[i] for i in range(len(timePoints)-1)]
    # initialize the return value with the first evaluation
    timeseries = [(0, initialValue)]

    #iterate over the relative time points and absolute time points
    for i in range(len(relativeTimePoints)):
        # solve the ODE-System with initial value = last evaluation and stopping time the current distance to the next
        # time step
        y = getSolutionAtTimeT(initialValue=initialValue, stoppingTime=relativeTimePoints[i], dt=dt,
                               parameters=parameters, method=method)
        initialValue = y.y # update initial value for next computation step
        # save the result and the correct time point
        timeseries += [(timePoints[i+1], list(y.y))]
    return timeseries

def computeError(timeseries, parameters, dt, method):
    """
    This method computes for a set of parameters for the ODE-system a solution and its least square error
    with respect to the time series
    :param timeseries: list of tuples [(t_1,p_1), ... , (t_n, p_n)]
    :param parameters: parameters for the ODE-System
    :param dt: size of timesteps
    :return: least square error of the solution and the time series
    """
    timePoints = [t for (t,p) in timeseries]
    points = [p for (t,p) in timeseries]
    computedTimeseries= getSolutionAtTimepoints(initialValue=points[0], timePoints=timePoints, dt=dt,
                                                parameters=parameters, method=method)
    computedPoints = [p for (t,p) in computedTimeseries]
    error = 0
    for i in range(len(timeseries)):
        error += distance.euclidean(points[i], computedPoints[i])**2
    return error

def plotODEsolution(parameters, initialValue, stoppingTime, numberOfSteps, method):
    """
    Plots the solution of the ODE-system loaded
    :param parameters: parameters in list form
    :param initialValue: initial state of the solution
    :param stoppingTime:
    :return: -
    """
    dt = stoppingTime/numberOfSteps
    y = ode(f).set_integrator(method)
    y.set_initial_value(initialValue, 0.0).set_f_params(parameters) # set initial value at time = 0
    evaluationTimes = [0.0] # initialized
    solution = [initialValue] # save the first time step
    while y.successful() and y.t < stoppingTime:
        evaluationTimes += [y.t+dt]
        y.integrate(y.t+dt)
        solution += [list(y.y)]
        if y.successful() is False:
            print("Something went wrong during r.integrate()")
    plt.ion()
    plt.axis([0.0, stoppingTime, 0.0, 1.1])
    for i in range(len(initialValue)):
        componentOfSolution = [solution[j][i] for j in range(len(solution))] # extract i-th component of solution vector
        plt.plot(evaluationTimes, componentOfSolution, label='x'+str(i+1))

    plt.ylabel('x')
    plt.xlabel('time')
    plt.legend(loc=0)
    plt.title("Trajectory of the solutions of the ODE-system with initial state "+str(initialValue))
    plt.show(block=True)



