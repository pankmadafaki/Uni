import numpy as np
### a)
def legendre(n,x):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return ((((2*n) - 1)*x*legendre(n-1,x)) - ((n-1)*legendre(n-2,x)))/(float(n))

### b)
def deriv_of_legendre(n,x):
    return (((n + 1)*x*legendre(n,x)) - ((n + 1)*legendre(n+1,x))) / (1 - (x**2))

n = 3
const = 0.0000000000000001
guess_value = -0.9
x_old = guess_value
zeros = []
def zeros_of_legendre():
    global x_old
    global i
    global guess_value
    global zeros
    i = 0
    while True:
        x_new = float(x_old - (legendre(n,x_old) / deriv_of_legendre(n,x_old)))
        difference = x_new - x_old
        if abs(difference) < const:
            zeros.append(x_new)
            i = 0
            guess_value = round(guess_value + 0.1,1)
            if guess_value == 0.0:
                guess_value = 0.1
            x_old = float(guess_value)
        else:
            x_old = x_new
            i = i + 1
        if guess_value == 1.0:
            zeros = list(dict.fromkeys(zeros))
            print("Zeros of the Legendre polynomial of degree",str(n),": ", zeros)
            break

zeros_of_legendre()