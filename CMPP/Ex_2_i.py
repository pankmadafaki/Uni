import matplotlib.pyplot as plt
import numpy as np

EN = 10
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


def function(x):
    return (0.2)/(x + 0.5) + 0.002



xk, zk = x(EN, 1)
print(xk, zk)
zk_half = np.zeros(EN-2)
xk_half = np.zeros(EN-2)
for k in range(EN-2):
    zk_half[k] = (zk[k]+zk[k+1])/2
    xk_half[k] = ((1+nthroot(zk_half[k], 1))/(1-nthroot(zk_half[k], 1)))




def linear_interpolation(x, z, x_half, z_half):
    f_k_half = np.zeros(len(z_half))
    for k in range(len(z_half)):
        f_k_half[k] = function(x[k]) + ((function(x[k+1])-function(x[k]))/(x[k+1]-x[k]))*(x_half[k] - x[k])
    return f_k_half

plt.scatter(xk, function(xk))
plt.scatter(xk_half, linear_interpolation(xk, zk, xk_half, zk_half))
plt.show()