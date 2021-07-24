import matplotlib.pyplot as plt
import numpy as np


def function(x):
    return (0.2)/(x + 0.5) + 0.002

def sum(x, M):
    summation = 0
    for i in range(1, M-1):
        summation += x
    return summation

#---------------------------------------------zk xk- create---------------------------------------------------------------
EN = 20
def nthroot(f, n):
    if f >= 0:
        return f**(1/n)
    else:
        return -abs(f)**(1/n)

def z(x, m):
    return (pow(x-1, m))/(pow(x+1, m))

def x(N, m):
    zk = np.zeros(N-1)
    xk = np.zeros(N-1)
    for k in range(N-1):
        zk[k] = -1 + 2 * (k+1) / N
        xk[k] = ((1+nthroot(zk[k], m))/(1-nthroot(zk[k], m)))
    return xk, zk

def T_k(x, k):
    if k == 0:
        temp = np.empty(len(x)); temp.fill(1)
        return temp
    elif k == 1:
        return x
    else:
        return 2*x*T_k(x, k-1) - T_k(x, k - 2)

def T_k_x(x, k):
    if k == 0:
        return 1
    elif k == 1:
        return x
    else:
        return 2*x*T_k_x(x, k-1) - T_k_x(x, k - 2)


def d_i_2(x_dat, i):
    N = len(x_dat)
    return np.sum(function(x_dat)*T_k(x_dat, i))*(2/N)

def d_i_cos(x_dat, i):
    N = len(x_dat)
    return np.sum(np.cos(x_dat)*T_k(x_dat, i))*(2/N)

def d_i(x_dat, y_dat, i):
    N = len(x_dat)
    return np.sum(y_dat*T_k(x_dat, i))*(2/N)


xk, zk = x(EN, 1)
#print(xk, zk)
zk_half = np.zeros(EN-1)
xk_half = np.zeros(EN-1)
for k in range(EN-2):
    zk_half[k] = (zk[k]+zk[k+1])/2
    xk_half[k] = ((1+nthroot(zk_half[k], 1))/(1-nthroot(zk_half[k], 1)))

fk = function(xk)



plt.scatter(xk, fk)
plt.scatter(xk_half, function(xk_half))
plt.show()


x_data = np.linspace(-1, 1, 20)
x_to_interpolate = np.linspace(-1, 1, 50)
y_data = np.cos(x_data*2*np.pi)
#y_data = x_data**3
#y_data = np.exp((-1)*x_data**2)
#y_to_interoplate = x_to_interpolate**3
y_to_interoplate = np.cos(x_to_interpolate*2*np.pi)
plt.scatter(x_data, y_data)
plt.scatter(x_to_interpolate, y_to_interoplate)
plt.show()


def chebyshev_interpolation(x_dat, y_dat, x_int_grid):
    m = int(len(x_dat)*0.9)

    y_interpolated = np.zeros(len(x_int_grid))
    for i in range(len(x_int_grid)):

        i_above, i_below = np.where(x_dat == x_dat[x_dat >= x_int_grid[i]].min())[0][0],\
                           np.where(x_dat == x_dat[x_dat <= x_int_grid[i]].max())[0][0]
        i_above = i_below
        summation = 0
        for j in range(1, m-1):
            summation += d_i(x_dat, y_dat, j)*T_k_x(x_int_grid[i], j)

        y_interpolated[i] = summation + d_i(x_dat, y_dat, 0)/2

    return y_interpolated


def clenshaw_recurrence(x_dat, y_dat, x_int_grid):
    m = int(len(x_dat)*0.9)
    y_interpolated = np.zeros(len(x_int_grid))
    print(m)
    print(len(x_dat), len(x_int_grid))

    for n in range(len(x_int_grid)):
        e = np.zeros(m+2)
        e[m] = 0
        e[m-1] = 0
        for j in range(m - 1, 0, -1):
            e[j] = 2*x_int_grid[n]*e[j + 1] - e[j + 2] + d_i(x_dat, y_dat, j)

        y_interpolated[n] = x_int_grid[n]*e[1] - e[2] + d_i(x_dat, y_dat, 0)/2

    return y_interpolated


chebyshev_interpolation(x_data, y_data, x_to_interpolate)

plt.plot(x_to_interpolate, chebyshev_interpolation(x_data, y_data, x_to_interpolate))
plt.scatter(x_data, y_data)
plt.show()
print("\n \n \n break \n")


clenshaw_recurrence(x_data, y_data, x_to_interpolate)

plt.plot(x_to_interpolate, clenshaw_recurrence(x_data, y_data, x_to_interpolate))
plt.scatter(x_data, y_data)
plt.show()
print("\n \n \n brealolk \n")

plt.plot(xk_half, clenshaw_recurrence(zk, fk, zk_half))
plt.scatter(xk, fk)
plt.show()




