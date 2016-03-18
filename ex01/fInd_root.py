# Bisection Method use in python3
import math


class fDivergence(Exception):
    pass


class kOverflow(Exception):
    pass


def diff(func, x):
    delta = x / 1000
    return (func(x + delta) - func(x - delta)) / (2.0 * delta)


def bolzano(func, left_number, right_number , f_error=1e-10 , k_max = 100):
    # define
    f = func
    l = left_number
    r = right_number
    k = 0

    if f(l) * f(r) > 0:
        raise fDivergence

    while abs(f(l) * f(r)) > f_error**2:
        t = (l + r) / 2
        if f(l) * f(t) > 0:
            l = t
        else:
            r = t
        k += 1
        if k > k_max:
            raise kOverflow
    if abs(f(l)) < f_error:
        if abs(f(r)) < f_error:
            root = (l + r) / 2
        else:
            root = l
    else:
        root = r
    return {'root': root, 'f_error': f(root),'x_error': abs(l-r),'k':k,'bolzano': "bolazno"}


def picard(func, guess_root, f_error=1e-10, k_max=100):
    f = func
    k = 1
    x1 = guess_root
    x2 = f(x1)

    # try solving divergence problem
    # print("abs((f(x2)-f(x1))/(x2-x1)) = ",abs((f(x2)-f(x1))/(x2-x1)))
    # if abs((f(x2) - f(x1)) / (x2 - x1)) > 1:
    #     raise fDivergence

    while abs(x2-x1) > f_error:
        (x1,x2) = (x2,f(x2))
        k += 1
        if k>k_max:
            raise kOverflow
    return {'root' : x2, 'f_error': f(x2) - x2,'x_eroor' : abs(x2-x1),'k':k, 'method': "picard"}


def newton(func, guess_root, f_error = 1e-10, k_max=100):
    f = func
    x1 = guess_root
    x2 = x1 - f(x1) / diff(func, x1)
    k = 1
    while abs(f(x2)) > f_error:
        x1 = x2
        x2 -= f(x2)/ diff(func,x2)
        k += 1
        if k > k_max:
            raise kOverflow
    return {'root' : x1,'f_error': f(x1),'x_error':abs(x2-x1),'k': k, 'method': "Newton"}


def newton_downhill(func, guess_root, f_error=1e-10, k_max=100):
    class qOverflow(Exception):
        pass
    f = func
    x1 = guess_root
    x2 = x1-f(x1)/diff(func,x1)
    k = 1
    while abs(f(x2)) > f_error:
        q = 1
        temp = x2
        while abs(f(x1)) <= abs(f(x2)):
            x2 = x1 - q * f(x1)/diff(f, x1)
            q /= 2
            if q == 2**-10:
                raise qOverflow
        x1 = temp
        k += 1
        if k > k_max:
            raise kOverflow
    return {'root': x2, 'f_error': f(x2), 'x_error': abs(x2 - x1), 'k': k, 'method':"Newton Downhill"}


def post_acceleration(func, guess_root, f_error=1e-10, k_max=100):
    f = func
    k = 1
    x1 = guess_root
    x2 = f(x1)

    while abs(x2 - x1) > f_error:
        L = diff(f,x2)
        (x1, x2) = (x2, (f(x2)-L*x2)/(1-L))
        k += 1
        if k > k_max:
            raise kOverflow
    return {'root': x2, 'f_error': f(x2) - x2, 'x_eroor': abs(x2 - x1), 'k': k, 'method': "post_acceleration"}


def aitken_acceleration(func, guess_root, f_error=1e-10, k_max=100):
    f = func
    k = 1
    x1 = guess_root
    temp1 = f(x1)
    temp2 = f(temp1)
    x2 = temp2 - (temp2-temp1)**2/(temp2-2*temp1 + x1)
    while abs(x2 - x1)> f_error:
        x1 = x2
        temp1 = f(x1)
        temp2 = f(temp1)
        x2 = temp2 - (temp2 - temp1) ** 2 / (temp2 - 2 * temp1 + x1)
        k += 1
        if k > k_max:
            raise kOverflow
    return {'root': x2, 'f_error': f(x2) - x2, 'x_eroor': abs(x2 - x1), 'k': k, 'method': "aitken_acceleration"}

def main():
    print(newton_downhill(func=lambda x: x * math.exp(x) - 1, guess_root=0.5))
    print(newton(func=lambda x: x*math.exp(x)-1, guess_root=0.5))
    print(picard(func=lambda x: math.exp(-x), guess_root=0.5))
    print(post_acceleration(func=lambda x:math.exp(-x),guess_root=0.5))
    print(aitken_acceleration(func=lambda x: math.exp(-x), guess_root=0.5))


if __name__ == '__main__':
    main()