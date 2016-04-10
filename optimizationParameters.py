__author__ = 'Robert Schwieger'

#------------
"""
These imports specify the model which is used
"""
from models import getBoundsA_normalized as getBounds
from models import convertToHumanReadableParametersA_normalized as convertToHumanReadableParameters
from models import continuousHomologA_normalized as continuousHomolog
from models import getInitialGuessA_normalized as getInitialGuessT
#------------
"""
Different algorithms for the optimization which work with bounds
"""
from scipy.optimize import brute
from scipy.optimize import basinhopping

# method='L-BFGS-B' # L-BFGS-B algorithm
# method='TNC' # truncated Newton (TNC) algorithm
method='SLSQP' #  Sequential Least SQuares Programming (SLSQP)


ftol = 1e-50 # Precision goal for the value of f in the stopping criterion.

"""
Parameters for the ODE-solver
"""
numberOfSteps = 10**4 # number of steps used for numerical solver of the ODEs
methodForODE = 'lsoda' # ODE solver

"""
Available data for fitting
"""
# time series used for optimization
timeseries = [(0,[0,0,0,0]), (6+2/3,[1,0,0,0]), (13+1/3,[1,1,0,0]), (20,[1,1,1,0]), (30,[1,1,1,1]),
              (40,[0,1,1,1]), (60,[0,0,1,1]), (90,[0,0,0,1]), (120,[0,0,0,0])]

"""
# Chose initial parameters. In case they were chosen randomly print them here also.
"""
# x0 = [3.34013773829254, 4.4054725759589, 1.6498446533954731, 4.803270099094165, 1.7879077880546044, 0.4590739709622539, 0.9911358947388088, 0.732498049647393, 0.245478685316471, 0.6498874437473751, 3.1017785174679435, 3.215085849863918, 9.763370261868282, 8.58305429706469]
"""
result = [ 3.38818778,  4.31806483,  1.33354283,  4.79752918,  1.42254256,
        0.49992379,  0.50004566,  0.50010122,  0.2303437 ,  0.48306036,
        3.10518132,  3.22540971,  9.76633977,  8.58238796]
min = 8.0001719989097815
"""
min_k = 1
max_k = 5
min_theta = 0
max_theta = 1
min_d = 0
max_d = 10

def getInitialGuess():
    return getInitialGuessT(min_k=min_k, max_k=max_k, min_theta=min_theta, max_theta=max_theta, min_d=min_d, max_d=max_d)

x0 = getInitialGuess()
numberOfOptimizations = 20
"""
k, theta, d = convertToHumanReadableParameters(x0)
print("initial guess: ")
print("Internal format: "+str(x0))
print("In human readable format")
print("k_initial: "+str(k))
print("theta_initial: "+str(theta))
print("d_initial: "+str(d))
"""