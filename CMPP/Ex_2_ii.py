import matplotlib.pyplot as plt
import numpy as np


def function(x):
    return (0.2)/(x + 0.5) + 0.002


p_data = np.array ( [ \
   [ -1.0, 1.00 ], \
   [ -0.8, 0.64 ], \
   [ -0.6, 0.36 ], \
   [ -0.4, 0.16 ], \
   [ -0.2, 0.04 ], \
   [  0.0, 0.00 ], \
   [  0.2, 0.04 ], \
   [  0.20001, 0.05 ], \
   [  0.4, 0.16 ], \
   [  0.6, 0.36 ], \
   [  0.8, 0.64 ], \
   [  1.0, 1.00 ] ] )

range_p_min = min(p_data[:,0])
range_p_max = max(p_data[:,0])
print(range_p_min, range_p_max)


#---------------------------------------------zk xk- create---------------------------------------------------------------
EN = 50
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



xk, zk = x(EN, 1)
print(xk, zk)
zk_half = np.zeros(EN-2)
xk_half = np.zeros(EN-2)
for k in range(EN-2):
    zk_half[k] = (zk[k]+zk[k+1])/2
    xk_half[k] = ((1+nthroot(zk_half[k], 1))/(1-nthroot(zk_half[k], 1)))

fk = function(xk)

#---------------------------------------------zk - create---------------------------------------------------------------



#----------------------------------------Tri-Diagonal_solver-----------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------


def A(x, x_interp, i):
    return (x[i + 1] - x_interp) / (x[i + 1] - x[i])


def B(x, x_interp, i):
    return (x_interp - x[i]) / (x[i + 1] - x[i])


def C(x, x_interp, i):
    return (1 / 6) * (A(x, x_interp, i) ** 3 - A(x, x_interp, i)) * (x[i + 1] - x[i]) ** 2


def D(x, x_interp, i):
    return (1 / 6) * (B(x, x_interp, i) ** 3 - B(x, x_interp, i)) * (x[i + 1] - x[i]) ** 2

def tridiagonal_solver(a, b, c, r):
    N = len(b)

    n = N - 1
    M = np.zeros((N, N))
    if len(a) == len(b) - 1 and len(c) == len(b) - 1:
        print("great  success")
    else:
        print("fail")
    beta = np.zeros(N)
    rho = np.zeros(N)
    z = np.zeros(N)
    beta[0], rho[0] = b[0], r[0]

    for i in range(1, N, 1):
        beta[i] = b[i] - (a[i-1]/beta[i-1])*c[i-1]
        rho[i] = r[i] - (a[i-1]/beta[i-1])*rho[i-1]


    z[N-1] = rho[N-1]/beta[N-1]
    for j in range(1, N, 1):
        z[N-1-j] = (rho[N-1-j]-c[N-1-j]*z[N-1-j+1])/(beta[N-1-j])


    return z

test_a = np.empty(9); test_a.fill(1)
test_b = np.empty(10); test_b.fill(2)
test_c = np.empty(9); test_c.fill(3)
test_r = np.empty(10); test_r.fill(4)

def generator(x, y):
    l = len(x)
    a = np.zeros(l - 3)
    b = np.zeros(l - 2)
    c = np.zeros(l - 3)
    r = np.zeros(l - 2)
    for i in range(1, l-1, 1):
        r[i-1] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1])
        b[i-1] = (x[i+1] - x[i-1])/3

    for i in range(1, l-2, 1):
        a[i - 1] = (x[i] - x[i - 1]) / 6
        c[i - 1] = (x[i + 1] - x[i]) / 6

    return a, b, c, r


def cubic_spline_interpolation(x_dat, y_dat, x_int_grid, z_dat):
    y_sec = np.zeros(len(x_dat))
    y_sec[1:len(x_dat)-1] = z_dat
    y_sec[0] = y_sec[1]
    y_sec[len(x_dat)-1] = y_sec[len(x_dat)-2]
    y_interpolated = np.zeros(len(x_int_grid))
    for i in range(len(x_int_grid)-1):
        print(i)
        i_above, i_below = np.where(x_dat == x_dat[x_dat >= x_int_grid[i]].min())[0][0],\
                           np.where(x_dat == x_dat[x_dat <= x_int_grid[i]].max())[0][0]
        i_above = i_below
        y_interpolated[i] = A(x_dat, x_int_grid[i], i_above)*y_dat[i_above]\
                            + B(x_dat, x_int_grid[i], i_above)*y_dat[i_above+1] \
                            + C(x_dat, x_int_grid[i], i_above)*y_sec[i_above] \
                            + D(x_dat, x_int_grid[i], i_above)*y_sec[i_above+1]
        print(y_interpolated[i])


    return y_interpolated



#x_data = np.linspace(0, 4*np.pi, 10)
x_data = p_data[:,0]
x_to_interpolate = np.linspace(range_p_min, range_p_max, 100)
print(np.shape(x_to_interpolate))
y_data = p_data[:,1]


a, b, c, r = generator(x_data, y_data)
z = tridiagonal_solver(a, b, c, r)

plt.plot(x_to_interpolate, cubic_spline_interpolation(x_data, y_data, x_to_interpolate, z))
#plt.plot(x_to_interpolate, np.tan(x_to_interpolate), color = 'r')
plt.scatter(x_data, y_data, edgecolors='black', color = 'r')
plt.show()

ak, bk, ck, rk = generator(xk, fk)
z_k = tridiagonal_solver(ak, bk, ck, rk)

plt.plot(xk_half, cubic_spline_interpolation(xk, fk, xk_half, z_k))
plt.scatter(xk, fk, edgecolors='black', color = 'r')
plt.show()



#print(np.argmin(x_dat[x_dat >= x_int_grid[i]]))
#        my_list = x_dat[x_dat >= x_int_grid[i]]
 #       val_above, idx_above = min((val_above, idx_above) for (idx_above, val_above) in enumerate(my_list))
 #       print(val_above, idx_above)