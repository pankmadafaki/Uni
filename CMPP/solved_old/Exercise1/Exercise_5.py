import numpy as np
x0 = 0
xn = 10
points = 100
x_i_array = np.around(np.linspace(x0, xn, points), 3)
interval = []

def function(y):
    return np.exp(y)

def check_function(x_i_array, xn, x):
    for i in range(len(x_i_array)):
        xi = round(x_i_array[i], 2)
        xi1 = round(x_i_array[i+1], 2)
        if i == 0 or i == points-1:
            if x < xi or x > xn:
                print("Error, bad idea to extrapolate! x is outside the interval.")
                break
        if x >= xi and x < xi1:
            interval.append(xi)
            interval.append(xi1)
            interval.append(x)
            print("x = ", x)
            print("Found the interval! The values are ", xi, "and" , xi1)
            break


def interpolation_formula(x_i_array, xn, x):
    check_function(x_i_array, xn, x)
    if len(interval) == 3:
        x_i, x_i1, x_value = interval[0], interval[1], interval[2]
        
        return result
    else:
        print("There was an error in the check function. Try again with different x value!")

my_value = interpolation_formula(x_i_array, xn, 0.3)
print(my_value)
