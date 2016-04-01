# solving linear_equation use in python
import numpy as np


class wrongShape(Exception):
    pass


class kOverflow(Exception):
    pass


def switchmatrix(flag="row", matrix=np.matrix(np.arange(12).reshape(3, 4)), m=1, n=0):
    if flag == "cloumn":
        (matrix[:, m], matrix[:, n]) = (matrix[:, n].copy(), matrix[:, m].copy())
    if flag == "row":
        (matrix[m, :], matrix[n, :]) = (matrix[n, :].copy(), matrix[m, :].copy())


def downtri(A, B):
    shape = A.shape[0]
    X = np.matrix(np.zeros((shape, 1)), dtype=np.double)
    for i in range(shape):
        X[i] = (B[i] - np.dot(A[i, :], X)) / A[i, i]
    return X


def uptri(A, B):
    shape = A.shape[0]
    X = np.matrix(np.zeros((shape, 1)), dtype=np.double)
    for i in range(shape - 1, -1, -1):
        X[i] = (B[i] - np.dot(A[i, :], X)) / A[i, i]
    return X


def gausse_limination(A, B):
    matrix = np.column_stack((A, B))
    if A.shape[0] - A.shape[1] != 0 or B.shape[1] - 1 != 0:
        raise wrongShape
    shape = A.shape[0]
    for i in range(shape):
        jmaxindex = int(i + np.argmax(matrix[i:, i]))
        switchmatrix(matrix=matrix, m=i, n=jmaxindex)
        matrix[i, :] /= matrix[i, i]
        for j in range(i + 1, shape):
            matrix[j, :] -= matrix[j, i] * matrix[i, :]
            j += 1
        i += 1
    return uptri(matrix[:, :-1], matrix[:, -1])


def lu(matrix=np.matrix(np.double(np.arange(3, 28).reshape(5, 5)))):
    if matrix.shape[0] - matrix.shape[1] != 0:
        raise wrongShape
    shape = matrix.shape[0]
    L = np.matrix(np.identity(shape))
    U = np.matrix(np.zeros([shape, shape]))
    U[0, :] = matrix[0, :]
    L[:, 0] = matrix[:, 0] / U[0, 0]
    for i in range(1, shape):
        for j in range(i, shape):
            U[i, j] = matrix[i, i] - np.dot(L[i, :], U[:, j])
        for j in range(i + 1, shape):
            L[j, i] = (matrix[j, i] - np.dot(L[j, :], U[:, i])) / U[i, i]
    return L, U


def doolittle(A, B):
    L, U = lu(A)
    return uptri(U, downtri(L, B))


def jacobi(A, B, guess_root, f_error=1e-3, k_max=100):
    if A.shape[0] - A.shape[1] != 0:
        raise wrongShape
    shape = A.shape[0]
    D = np.matrix(np.diag(np.diag(A)))
    G = np.identity(shape) - D.I * A
    g = D.I * B
    X = G * np.matrix(guess_root) + g
    k = 1
    while np.linalg.norm(A * X - B) > f_error:
        X0 = X
        X = G * X + g
        k += 1
        if k > k_max:
            raise kOverflow
    return {'root': X, 'f_error': np.linalg.norm(A * X - B), 'x_eroor': np.linalg.norm(X - X0), 'k': k,
            'method': "jacobi"}


def gauss_seidel(A, B, guess_root, f_error=1e-3, k_max=100):
    if A.shape[0] - A.shape[1] != 0:
        raise wrongShape
    shape = A.shape[0]
    X = np.matrix(guess_root)
    k = 0
    while np.linalg.norm(A * X - B) > f_error:
        for i in range(shape):
            X0 = X.copy()
            X_copy = X.copy()
            X_copy[i, 0] = 0
            X[i, 0] = -(A[i, :] * X_copy - B[i, 0])[0, 0] / A[i, i]
        k += 1
        if k > k_max:
            raise kOverflow
    return {'root': X, 'f_error': np.linalg.norm(A * X - B), 'x_eroor': np.linalg.norm(X - X0), 'k': k,
            'method': "gauss_seidel"}


def successive_relaxation(A, B, guess_root, omega=1.4, f_error=1e-9, k_max=10000):
    if A.shape[0] - A.shape[1] != 0:
        raise wrongShape
    shape = A.shape[0]
    X = np.matrix(guess_root)
    k = 0
    while np.linalg.norm(A * X - B) > f_error:
        for i in range(shape):
            X0 = X.copy()
            X_copy = X.copy()
            x_i = (1 - omega) * X_copy[i, 0]
            X_copy[i, 0] = 0
            X[i, 0] = x_i - omega * (A[i, :] * X_copy - B[i, 0])[0, 0] / A[i, i]
        k += 1
        if k > k_max:
            raise kOverflow
    return {'root': X, 'f_error': np.linalg.norm(A * X - B), 'x_eroor': np.linalg.norm(X - X0), 'k': k,
            'method': "successive_relaxation"}


def main():
    A = np.matrix(np.double([[31, -13, 0, 0, 0, -10, 0, 0, 0],
                             [-13, 35, -9, 0, -11, 0, 0, 0, 0],
                             [0, -9, 31, -10, 0, 0, 0, 0, 0],
                             [0, 0, -10, 79, -30, 0, 0, 0, -9],
                             [0, 0, 0, -30, 57, -7, 0, -5, 0],
                             [0, 0, 0, 0, -7, 47, -30, 0, 0],
                             [0, 0, 0, 0, 0, -30, 41, 0, 0],
                             [0, 0, 0, 0, -5, 0, 0, 27, -2],
                             [0, 0, 0, -9, 0, 0, 0, -2, 29]
                             ]))
    B = np.matrix(np.double([[-15],
                             [27],
                             [-23],
                             [0],
                             [-20],
                             [12],
                             [-7],
                             [7],
                             [10]
                             ]))
    X = np.matrix(np.double([[0],
                             [0],
                             [0],
                             [0],
                             [0],
                             [0],
                             [0],
                             [0],
                             [0]
                             ]))

    Y = successive_relaxation(A, B, X, omega=1.5, f_error=1e-3)


if __name__ == '__main__':
    main()