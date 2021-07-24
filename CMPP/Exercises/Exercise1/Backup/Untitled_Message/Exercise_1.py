import numpy as np
a = 0
b = 10
n = 1000
x = np.linspace(a,b,n+1)
h = float((b-a)/n)

def f_1(x):
    return np.sin(x)*np.e**(-(np.sqrt(x)))
def f_2(x):
    return np.tanh(x)

### Trapezoidal rule
sum_number_f1 = 0
sum_number_f2 = 0
for i in x:
    if i == x[0] or i == x[-1]:
        f1_i = f_1(i)
        f2_i = f_2(i)
    else:
        f1_i = 2*f_1(i)
        f2_i = 2*f_2(i)
    sum_number_f1 = sum_number_f1 + f1_i
    sum_number_f2 = sum_number_f2 + f2_i

result_I1 = (h/2) * sum_number_f1
result_I2 = (h/2) * sum_number_f2
print('Result of first function with trapezoid rule: ', result_I1)
print('Result of second function with trapezoid rule: ', result_I2)
### For the first function, at around 630 points we get the first 4 digits correct (0.503)
### For the second function, at around 100 points we get the first 4 digits correct (9.308)

### Simpsons one third rule
#print('x: ', x)
sum_number2_f1 = 0
sum_number2_f2 = 0
for i in range(len(x)):
    list_elem = x[i]
    if list_elem == x[0] or list_elem == x[-1]:
        f1_i = f_1(list_elem)
        f2_i = f_2(list_elem)
    elif i % 2 == 0:
        f1_i = 2 * f_1(list_elem)
        f2_i = 2 * f_2(list_elem)
    else:
        f1_i = 4 * f_1(list_elem)
        f2_i = 4 * f_2(list_elem)
    sum_number2_f1 = sum_number2_f1 + f1_i
    sum_number2_f2 = sum_number2_f2 + f2_i

result2_I1 = (h/3) * sum_number2_f1
result2_I2 = (h/3) * sum_number2_f2
print('Result of first function with Simpsons 1/3 rule: ', result2_I1)
print('Result of second function with Simpsons 1/3 rule: ', result2_I2)
### For the first function, at around 150 points we get the first 4 digits correct (0.503)
### For the second function, at around 50 points we get the first 4 digits correct (9.308)