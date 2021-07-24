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


n = 20
const = float(1.2 * 10 ** -16)
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
            guess_value = float(np.cos((np.pi * (i - 0.25)) / (n + 0.5)))
            # print("guess value: ", guess_value)
            x_old = guess_value
        else:
            x_old = x_new
        if len(zeros) == n:
            for j in range(n):
                w_j = 2 / ((1 - (zeros[j] ** 2)) * ((deriv_of_legendre(n, zeros[j])) ** 2))
                weights.append(w_j)
            # print("Zeros of the Legendre polynomial of degree",str(n),": ", zeros)
            # print("Weights of the zeros for Legendre polynomial of degree", str(n), ": ", weights)
            break


zeros_and_weigths_of_legendre()
print("Zeros", zeros)
print("Weights", weights)

### b)
x_min = -1
x_max = 1
a_magnitude = 0
b_magnitude = 1
a_factor = 5
b_factor = 5
a = float(a_factor * 10 ** a_magnitude)
print(a)
b = float(b_factor * 10 ** b_magnitude)
print(b)
x_points = np.linspace(x_min, x_max, 100)
z_points = np.linspace(a, b, 100)
A = (-np.log(a * b)) / (np.log(b / a))
B = 2 / (np.log(b / a))
print(A, B)

def log_transformation():
    return A + (B * np.log(z_points))

def inverse_log_transformation():
    return np.exp((x_points - A) / (B))

### c
result = 0
jacobian = (np.exp((zeros - A) / (B))) / (B)

def inverse_log_transformation2(x_i):
    return np.exp((x_i - A) / (B))


def function(x_i):
    return np.exp(inverse_log_transformation2(x_i)) - inverse_log_transformation2(x_i)
    #return np.sin(inverse_log_transformation2(x_i)) * np.e ** (-(np.sqrt(inverse_log_transformation2(x_i))))
    # return np.tanh(inverse_log_transformation(x_i))

def gauss_legendre_integration():
    global result
    for k in range(n):
        x_i = zeros[k]
        result = result + (jacobian[k] * weights[k] * function(x_i))
    print("Final result: ", result)


gauss_legendre_integration()
