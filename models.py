__author__ = 'Robert Schwieger'
from random import uniform
import math

def hillFunction(x, k, theta):
    # if the x-value is minimally negative it is set back to zero. However when the x-value for some reason falls
    # below neg_tol an Exception is thrown
    neg_tol = -10**(-2)
    if x == 0:
        return 0
    if x < 0:
        if x > neg_tol:
            print("x = "+str(x)+" < 0 in Hill function. Continue With x=0")
            return 0
        else:
            print("x = "+str(x)+" is negative.")
            raise Exception
    try:
        result = math.pow(x,k)/(math.pow(x,k)+math.pow(theta,k))
    except ValueError:
        print("x: "+str(x))
        print("k: "+str(k))
        print("theta: "+str(theta))
        raise ValueError
    return result

"""
normalized A model
"""
def continuousHomologA_normalized(x, k, theta):
    """
    odefy transformation
    :param x: list of real numbers in the interval [0,1]
    :param k: matrix of positive real numbers
    :param theta: matrix of real numbers in the interval (0,1)
    :return: list of real numbers in the interval [0,1]
    """
    try:
        y = [0,0,0,0]
        y[0] = 1-hillFunction(x[2], k[0][2], theta[0][2]) # P <- 1-A
        y[1] = hillFunction(x[0], k[1][0], theta[1][0]) # B <- P
        y[2] = hillFunction(x[1], k[2][1], theta[2][1]) # A <- B
        y[3] = hillFunction(x[2], k[3][2], theta[3][2])*hillFunction(x[0], k[3][0], theta[3][0]) # Ap <- A*P
    except ValueError:
        print("Value Error in continuousHomolog")
        raise ValueError
    return y

def convertToHumanReadableParametersA_normalized(parameters):
    """
    The parameters are required to be a list. To keep the code more human readable they are converted into the
    usual data structure used in the text.
    :param parameters: list of parameters
    :return: (k, theta, d)
    """
    k = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]
    k[0][2] = parameters[0]
    k[1][0] = parameters[1]
    k[2][1] = parameters[2]
    k[3][0] = parameters[3]
    k[3][2] = parameters[4]
    theta = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]
    theta[0][2] = parameters[5]
    theta[1][0] = parameters[6]
    theta[2][1] = parameters[7]
    theta[3][0] = parameters[8]
    theta[3][2] = parameters[9]
    d = [parameters[10], parameters[11], parameters[12], parameters[13]]
    return k, theta, d

def getInitialGuessA_normalized(min_k, max_k, min_theta, max_theta, min_d, max_d):
    return [uniform(1,5), uniform(1,5), uniform(1,5), uniform(1,5) ,uniform(1,5), # k
              uniform(0,1), uniform(0,1), uniform(0,1), uniform(0,1), uniform(0,1), # theta
              uniform(0,10),uniform(0,10),uniform(0,10),uniform(0,10)] # d


def getBoundsA_normalized():
    """
    Returns the bounds of the parameter for the optimization
    :return:
    """
    eps = 10**(-3)
    return [(1+eps, None), (1+eps, None), (1+eps, None), (1+eps, None), (1+eps, None), # Hill coefficients
        (eps,1-eps), (eps,1-eps), (eps,1-eps), (eps,1-eps), (eps,1-eps), # Thresholds
        (eps,None), (eps,None), (eps,None), (eps,None)]

"""
normalized B model
"""

def continuousHomologB_normalized(x, k, theta):
    """
    odefy transformation
    :param x: list of real numbers in the interval [0,1]
    :param k: matrix of positive real numbers
    :param theta: matrix of real numbers in the interval (0,1)
    :return: list of real numbers in the interval [0,1]
    """
    try:
        y = [0,0,0,0]
        y[0] = 1 # P <- 1 # ToDo: Noch mal anschauen. Muss man hier noch einen Parameter hinzufügen?
        y[1] = (1-hillFunction(x[2], k[1][2], theta[1][2]))*(1-hillFunction(x[1], k[1][1], theta[1][1])) # B <- (1-A)*(1-B)
        y[2] = hillFunction(x[1], k[2][1], theta[2][1]) # A <- B
        y[3] = hillFunction(x[2], k[3][2], theta[3][2])*hillFunction(x[0], k[3][0], theta[3][0]) # Ap <- A*P
    except ValueError:
        print("Value Error in continuousHomolog")
        raise ValueError
    return y

def convertToHumanReadableParametersB_normalized(parameters):
    """
    The parameters are required to be a list. To keep the code more human readable they are converted into the
    usual data structure used in the text.
    :param parameters: list of parameters
    :return: (k, theta, d)
    """
    k = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]
    k[1][1] = parameters[0]
    k[1][2] = parameters[1]
    k[2][1] = parameters[2]
    k[3][0] = parameters[3]
    k[3][2] = parameters[4]
    theta = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]
    theta[1][1] = parameters[5]
    theta[1][2] = parameters[6]
    theta[2][1] = parameters[7]
    theta[3][0] = parameters[8]
    theta[3][2] = parameters[9]
    d = [parameters[10], parameters[11], parameters[12], parameters[13]]
    return k, theta, d

getInitialGuessB_normalized = getInitialGuessA_normalized
getBoundsB_normalized = getBoundsA_normalized

"""
normalized AB model
"""

def continuousHomologAB_normalized(x, k, theta):
    """
    odefy transformation
    :param x: list of real numbers in the interval [0,1]
    :param k: matrix of positive real numbers
    :param theta: matrix of real numbers in the interval (0,1)
    :return: list of real numbers in the interval [0,1]
    """
    try:
        y = [0,0,0,0]
        y[0] = 1 # P <- 1
        y[1] = hillFunction(x[0], k[1][0], theta[1][0])*(1-hillFunction(x[3], k[1][3], theta[1][3])) # B <- P*(1-Ap)
        y[2] = hillFunction(x[1], k[2][1], theta[2][1]) # A <- B
        y[3] = hillFunction(x[2], k[3][2], theta[3][2])*hillFunction(x[0], k[3][0], theta[3][0]) # Ap <- A*P
    except ValueError:
        print("Value Error in continuousHomolog")
        raise ValueError
    return y