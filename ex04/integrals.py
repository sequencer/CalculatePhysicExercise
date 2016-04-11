import numpy as np
import matplotlib.pyplot as mpl
import math


def simpsonIntegral(function, xmin, xmax, N):
    dx = (xmax - xmin) / N
    X = []
    Y = []
    for i in range(N + 1):
        X.append(xmin + i * dx)
        Y.append(function(xmin + i * dx))
    s = Y[0] + Y[-1]
    for i in range(1, N - 1):
        if i % 2 != 0:
            s += 2 * function(X[i])
        else:
            s += 4 * function(X[i])
    s *= dx / 3
    return s


def trapezoidIntegral(function, xmin, xmax, N):
    dx = (xmax - xmin) / N
    X = []
    Y = []
    for i in range(N + 1):
        X.append(xmin + i * dx)
        Y.append(function(xmin + i * dx))
    s = Y[0] + Y[-1]
    for i in range(1, N - 1):
        s += 2 * function(X[i])
    s *= dx / 2
    return s


if __name__ == '__main__':
    function = lambda x: math.sin(x)
    xmin = 1
    xmax = 5
    N = 40
    real = np.cos(1) - np.cos(5)
    s = trapezoidIntegral(function, xmin, xmax, N)
    s = simpsonIntegral(function, xmin, xmax, N)
    trapezoiderror = []
    simpsoneroor = []
    X = list(range(50,100, 1))
    for i in X:
        trapezoiderror.append(trapezoidIntegral(function, xmin, xmax, i) - real)
        simpsoneroor.append(simpsonIntegral(function, xmin, xmax, i) - real)
        pass
    mpl.title('Error-N')
    mpl.xlabel('N')
    mpl.ylabel('Error->log(Error)')
    mpl.text(90,-2.6,'trapezoiderror', color='green')
    mpl.text(90,-2.7,'simpsoneroor',color = 'red')
    mpl.plot(X, np.log(trapezoiderror), color='green')
    mpl.plot(X, np.log(simpsoneroor), color='red')
    mpl.show()
