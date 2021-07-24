import numpy as np
import matplotlib.pyplot as plt
from scipy.special import eval_legendre



def legendre(n, x):
    return eval_legendre(n, x)


def deriv_of_legendre(n, x):
    return (((n + 1) * x * legendre(n, x)) - ((n + 1) * legendre(n + 1, x))) / (1 - (x ** 2))





def function(x_i):
    return np.sqrt(1-(x_i**2))


N = 100

def zeros_and_weigths_of_legendre(n):
    const = float(1.2 * 10 ** - 12)
    guess_value = float(np.cos((np.pi * (1 - 0.25)) / (n + 0.5)))
    x_old = guess_value
    zeros = []
    weights = []
    i = 1
    while True:
        x_new = float(x_old - (legendre(n, x_old) / deriv_of_legendre(n, x_old)))
        difference = x_new - x_old
        # print(difference)
        if abs(difference) < const:
            zeros.append(x_new)
            i = i + 1
            #print("New grid point ", i)
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
    result = 0
    for k in range(n):
        x_i = zeros[k]
        result = result + (weights[k] * function(x_i))
    print("Final result: ", result)
    return result



error = np.zeros(N)

for k in range(1, N, 1):
    error[k] = zeros_and_weigths_of_legendre(k) - np.pi/2
    print(k)

plt.plot(error, label= f"Final Error {error[N-1]}")
plt.legend()
plt.show()


