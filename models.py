__author__ = 'Robert Schwieger'
from random import uniform
import math

def hillFunction(x, k, theta):
    if x < 0:
        print("x<0 in Hill function")
        raise ValueError
    try:
        result = math.pow(x,k)/(math.pow(x,k)+math.pow(theta,k))
    except ValueError:
        print("x: "+str(x))
        print("k: "+str(k))
        print("theta: "+str(theta))
        print(math.pow(x,k))
        print(math.pow(theta,k))
        raise ValueError
    return result

def continuousHomologA_normalized(x, k, theta):
    """
    odefy transformation
    :param x: list of real numbers in the interval [0,1]
    :param k: matrix of positive real numbers
    :param theta: matrix of real numbers in the interval (0,1)
    :return: list of real numbers in the interval [0,1]
    """
    y = [0,0,0,0]
    y[0] = 1-hillFunction(x[2], k[0][2], theta[0][2]) # 1-A
    y[1] = hillFunction(x[0], k[1][0], theta[1][0]) # P
    y[2] = hillFunction(x[1], k[2][1], theta[2][1]) # B
    y[3] = hillFunction(x[2], k[3][2], theta[3][2])*hillFunction(x[0], k[3][0], theta[3][0])
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

def getInitialGuessA_normalized():
    return [uniform(1,5), uniform(1,5), uniform(1,5), uniform(1,5) ,uniform(1,5), # k
              uniform(0,1), uniform(0,1), uniform(0,1), uniform(0,1), uniform(0,1), # theta
              uniform(0,10),uniform(0,10),uniform(0,10),uniform(0,10)] # d

def constraintsSatisfiedA_normalized(parameters):
    """
    Checks if the constraints are satisfied
    :param parameters:
    :return: True = constraints satisfied, False = constraints are not satisfied
    """
    for i in range(5): # Hill coefficients
        if parameters[i]<=1:
            return False
    for i in range(5,10): # Thresholds
        if parameters[i]<=0 or parameters[i]>=1:
            return False
    for i in range(10,14): # D
        if parameters[i] <= 0:
            return False
    return True

def getBoundsA_normalized():
    """
    Returns the bounds of the parameter for the optimization
    :return:
    """
    eps = 10**(-3)
    return [(1+eps, None), (1+eps, None), (1+eps, None), (1+eps, None), (1+eps, None), # Hill coefficients
        (eps,1-eps), (eps,1-eps), (eps,1-eps), (eps,1-eps), (eps,1-eps), # Thresholds
        (eps,None), (eps,None), (eps,None), (eps,None)]