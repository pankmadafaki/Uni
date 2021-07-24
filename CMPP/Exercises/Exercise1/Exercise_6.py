import numpy as np

### a)
def chebyshev_polynomial(x,n):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return (2*x*(chebyshev_polynomial(x, n-1)) - chebyshev_polynomial(x, n - 2))

print(chebyshev_polynomial(3,2))

### b)
coeffictients = []
def coefficients():
    for i in range():
        res = 0
        for k in range(N-1):
            res = res + (funct(np.cos((np.pi * (k + 0.5)) / N)) * np.cos((np.pi * i * (k + 0.5)) / (N)) )
        d_i = res
        coeffictients.append(d_i)