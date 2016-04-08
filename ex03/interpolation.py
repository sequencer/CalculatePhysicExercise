import numpy as np
import matplotlib.pyplot as mpl
from linear_solving import downtri
from Function import linear_solving


class outOfBound(Exception):
    pass


def samplingFromFunction(function, xmin, xmax, N):
    X = np.linspace(xmin, xmax, N)
    Y = function(X)
    return X, Y


def larangeInterpolation(X, Y, X0):
    Y0 = []
    size = len(X)
    for x in X0:
        y = 0
        for i in range(size):
            li = 1
            for j in range(size):
                if j == i:
                    continue
                li *= (x - X[j]) / (X[i] - X[j])
            y += Y[i] * li
        Y0.append(y)
    return Y0


def newtonInterpolation(X, Y, X0):
    Y0 = []
    size = len(X)
    for x0 in X0:
        A = np.matrix(np.zeros((size - 1, size - 1)))
        for i in range(size - 1):
            for j in range(i + 1):
                A[i, j] = X[i + 1] - X[j]
                for k in range(j):
                    A[i, j] *= X[i + 1] - X[k]
        B = (np.matrix(Y[1:]) - Y[0]).transpose()
        M = np.matrix(np.zeros((1, size - 1)))

        for i in range(size - 1):
            M[0, i] = x0 - X[i]
            for j in range(i):
                M[0, i] *= x0 - X[j]
        print(M)
        Y0.append(np.double(M * downtri(A, B) + Y[0]))
    return Y0


def splineCurve(X, Y, X0):
    if X0[0] < X[0] or X0[-1] > X[-1]:
        raise outOfBound
    shape = len(X)
    A = 2 * np.identity(shape - 2)
    A = np.c_[np.zeros((shape - 2, 1)), A, np.zeros((shape - 2, 1))]
    mu = np.identity(shape - 2)
    la = np.identity(shape - 2)
    h = []
    C = np.zeros((shape - 2, 1), np.double)
    # 先做单位矩阵计算
    for i in range(shape - 1):
        h.append(X[i + 1] - X[i])
    for i in range(shape - 2):
        mu[i, i] = h[i+1] / (h[i] + h[i + 1])
        la[i, i] = 1 - mu[i, i]
        C[i] = 3 * (la[i, i] * (Y[i+1] - Y[i]) / h[i] + mu[i, i] * (Y[i + 2] - Y[i+1]) / h[i+1])
    # 添加行列后合并

    mu = np.c_[np.zeros((shape - 2, 2)), mu]
    la = np.c_[la, np.zeros((shape - 2, 2))]
    A += mu + la
    A = A[:,1:-1]
    A = np.matrix(A)
    M = linear_solving.gausse_limination(A, C)
    M = np.vstack((0,M,0))

    Y0 = []
    i = 0
    for x0 in X0:
        # 重复迭代让x0处于正确的范围内
        while X[i + 1] < x0:
            i += 1
            continue
        # 计算y0的值,如果越界使x+=1
        if X[i + 1] < x0:
            i +=  1
        # 计算y0
        y0 = 1 / h[i] ** 3 * (h[i] + 2 * (x0 - X[i])) * (x0 - X[i + 1]) ** 2 * Y[i] + 1 / h[i] ** 3 * (h[i] - 2 * (x0 - X[i + 1])) * (x0 - X[i]) ** 2 * Y[i + 1] + 1 / h[i] ** 2 * (x0 - X[i]) * (x0 - X[i + 1]) ** 2 * M[i, 0] + 1 / h[i] ** 2 * (x0 - X[i + 1]) * (x0 - X[i]) ** 2 * M[i + 1, 0]
        Y0.append(y0)
    return Y0



if __name__ == '__main__':
    function = lambda x: (1 + x ** 2) ** -1
    X, Y = samplingFromFunction(function, -5, 5, 15)
    X0 = np.linspace(-5, 5, 150)
    Y0 = splineCurve(X, Y, X0)
    mpl.plot(X0,Y0,'g')
    mpl.plot(X0,function(X0)-Y0,'b')
    mpl.show()
