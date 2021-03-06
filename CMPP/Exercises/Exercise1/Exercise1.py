import numpy as np
import matplotlib.pyplot as plt

### a)
def legendre(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return ((((2 * n) - 1) * x * legendre(n - 1, x)) - ((n - 1) * legendre(n - 2, x))) / (float(n))


def deriv_of_legendre(n, x):
    return (((n + 1) * x * legendre(n, x)) - ((n + 1) * legendre(n + 1, x))) / (1 - (x ** 2))


n = 50
const = float(1.2 * 10 ** - 12)
guess_value = float(np.cos((np.pi * (1 - 0.25)) / (n + 0.5)))
x_old = guess_value
zeros = []
weights = []


def zeros_and_weigths_of_legendre():
    global x_old
    global zeros
    global weights
    i = 1
    while True:
        x_new = float(x_old - (legendre(n, x_old) / deriv_of_legendre(n, x_old)))
        difference = x_new - x_old
        # print(difference)
        if abs(difference) < const:
            zeros.append(x_new)
            i = i + 1
            print("New grid point ", i)
            guess_value = float(np.cos((np.pi * (i - 0.25)) / (n + 0.5)))
            # print("guess value: ", guess_value)
            x_old = guess_value
        else:
            x_old = x_new
            print("Grid point not found, looking again ", i)
        if len(zeros) == n:
            for j in range(n):
                print(j)
                w_j = 2 / ((1 - (zeros[j] ** 2)) * ((deriv_of_legendre(n, zeros[j])) ** 2))
                weights.append(w_j)
            # print("Zeros of the Legendre polynomial of degree",str(n),": ", zeros)
            # print("Weights of the zeros for Legendre polynomial of degree", str(n), ": ", weights)
            break


zeros_and_weigths_of_legendre()
#print("Zeros", zeros)
#print("Weights", weights)

### c
result = 0

def function(x_i):
    return np.sqrt(1-(x_i**2))

def gauss_legendre_integration():
    global result
    for k in range(n):
        x_i = zeros[k]
        result = result + (weights[k] * function(x_i))
    print("Final result: ", result)

def IntegralResult():
    val = np.pi/2
    print("Integral result: ", val)

gauss_legendre_integration()
IntegralResult()