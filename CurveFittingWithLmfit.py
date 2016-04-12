__author__ = 'Robert Schwieger'

from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math
from random import uniform

"""
A model
"""
def hillFunction(x, k, theta):
    return math.pow(x,k)/(math.pow(x,k)+math.pow(theta,k))

def f(xs, t, ps):
    theta13 = ps['theta13'].value
    theta21 = ps['theta21'].value
    theta32 = ps['theta32'].value
    theta41 = ps['theta41'].value
    theta43 = ps['theta43'].value

    n13 = ps['n13'].value
    n21 = ps['n21'].value
    n32 = ps['n32'].value
    n41 = ps['n41'].value
    n43 = ps['n43'].value

    gamma1 = ps['gamma1'].value
    gamma2 = ps['gamma2'].value
    gamma3 = ps['gamma3'].value
    gamma4 = ps['gamma4'].value
    P, B, A, Ap = xs
    return [gamma1*(1-hillFunction(A, n13, theta13)-P), gamma2*(hillFunction(P, n21, theta21)-B),
            gamma3*(hillFunction(B, n32, theta32)-A), gamma4*(hillFunction(P, n41,theta41)*hillFunction(A, n43, theta43)-Ap)]

def g(t, x0, ps):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(f, x0, t, args=(ps,))
    return x

def residual(ps, ts, data):
    x0 = ps['P0'].value, ps['B0'].value, ps['A0'].value, ps['Ap0'].value
    model = g(ts, x0, ps)
    return (model - data).ravel()



def optimize(t, data):
    # set parameters incluing bounds
    params = Parameters()
    params.add('P0', value=float(data[0, 0]), min=0, max=1)
    params.add('B0', value=float(data[0, 1]), min=0, max=1)
    params.add('A0', value=float(data[0, 2]), min=0, max=1)
    params.add('Ap0', value=float(data[0, 3]), min=0, max=1)
    params.add('gamma1', value=uniform(0,10), min=0, max=10)
    params.add('gamma2', value=uniform(0,10), min=0, max=10)
    params.add('gamma3', value=uniform(0,10), min=0, max=10)
    params.add('gamma4', value=uniform(0,10), min=0, max=10)

    params.add('n13', value=uniform(1,10), min=1, max=10)
    params.add('n21', value=uniform(1,10), min=1, max=10)
    params.add('n32', value=uniform(1,10), min=1, max=10)
    params.add('n41', value=uniform(1,10), min=1, max=10)
    params.add('n43', value=uniform(1,10), min=1, max=10)

    params.add('theta13', value=uniform(0,1), min=0, max=1)
    params.add('theta21', value=uniform(0,1), min=0, max=1)
    params.add('theta32', value=uniform(0,1), min=0, max=1)
    params.add('theta41', value=uniform(0,1), min=0, max=1)
    params.add('theta43', value=uniform(0,1), min=0, max=1)

    # fit model and find predicted values
    result = minimize(residual, params, args=(t, data), method='leastsq') # leastsq "=" Levenberg-Marquardt
    return result

def optimizeNTimes(t, data, number):
    bestResult = None
    for i in range(number):
        print(str(i)+" of "+str(number))
        try:
            result = optimize(t, data)
        except:
            continue
        if bestResult is None or result.chisqr<bestResult.chisqr:
            bestResult = result
    return bestResult

def plotData(t, data, result):
    # plot data and fitted curves
    plt.axis([0.0, t[-1], -0.1, 1.1])
    plt.plot(t, data, 'o')
    plt.gca().set_color_cycle(None)
    # plt.plot(t, final, 'x', linewidth=2);

    t = np.linspace(0, 120, 100)
    data = g(t, data[0], result.params)
    plt.plot(t, data, '-', linewidth=2);
    plt.show()

def plotData2(t, data, result):
    numberOfComponents = 4
    f, axarr = plt.subplots(numberOfComponents, sharex=True)
    continuoust = np.linspace(0, 120, 100)
    continuousData = g(continuoust, data[0], result.params)

    for i in range(numberOfComponents):
        axarr[i].set_title('Component '+str(i+1))

        # plot data and fitted curves
        component = data[:, [i]] # print i-th component only
        axarr[i].axis([0.0, t[-1], -0.1, 1.1])
        axarr[i].plot(t, component, 'o')
        # plt.gca().set_color_cycle(None)


        component = continuousData[:, [i]] # print i-th component only
        axarr[i].plot(continuoust, component, '-', linewidth=2);

    plt.xlabel('time')
    plt.legend(loc=0)
    plt.show()

# time series data
t = [0, 6+2/3, 13+1/3, 20, 30, 40, 60, 90, 120]
data = np.array([[0, 0, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1],
        [0, 0, 0, 0]])

result = optimizeNTimes(t=t, data=data, number=1)
# display fitted statistics
report_fit(result)
#plot solution curve
print("residual error: "+str(result.residual))
print("chi-square: "+str((result.chisqr)))
print("chi: "+str((result.chisqr**0.5)))

plotData2(t, data, result)
