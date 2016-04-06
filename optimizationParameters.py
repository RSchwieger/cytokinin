__author__ = 'Robert Schwieger'

#------------
"""
These imports specify the model which is used
"""
from models import getBoundsA_normalized as getBounds
from models import convertToHumanReadableParametersA_normalized as convertToHumanReadableParameters
from models import continuousHomologA_normalized as continuousHomolog
from models import getInitialGuessA_normalized as getInitialGuess
#------------
"""
Different algorithms for the optimization which work with bounds
"""

# method='L-BFGS-B' # L-BFGS-B algorithm
# method='TNC' # truncated Newton (TNC) algorithm
method='SLSQP' #  Sequential Least SQuares Programming (SLSQP)

#-----------------


numberOfSteps = 10**4 # number of steps used for numerical solver of the ODEs
methodForODE = 'lsoda' # ODE solver

# time series used for optimization
timeseries = [(5,[0,0,0,0]), (10,[1,0,0,0]), (15,[1,1,0,0]), (20,[1,1,1,0]), (30,[1,1,1,1]),
              (40,[0,1,1,1]), (60,[0,0,1,1]), (90,[0,0,0,1]), (120,[0,0,0,0])]