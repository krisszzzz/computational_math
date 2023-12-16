import numpy as np
import matplotlib.pyplot as plt
from math import *

def f(x):
    if x == 0:
        return x
    else:
        return np.exp(-1 / x) / x**2 

step = 0.125
x = np.arange(0, 1 + step, step, dtype=np.float64)
f_val = [ f(x_val) for x_val in x ]
#print(x[8])

def trapezoidal(x, f):
    assert len(x) == len(f), 'Arrays must have the same length'
    I = 0
    for i in range(len(x) - 1):
        I = I + (x[i+1] - x[i]) / 2 * (f[i+1] + f[i])
    return I

def richardson(x, f, I, p):
    assert len(x) == len(f), 'Arrays must have the same length'
    return I(x, f) + (I(x, f) - I(x[::2], f[::2])) / (2**p - 1)

def simpson(x, f):
    assert len(x) == len(f), 'Arrays must have the same length'
    I = 0
    for i in range(len(x) // 2):
        I = I + (x[2*i+1] - x[2*i]) / 3 * (f[2*i] + 4*f[2*i+1] + f[2*i+2])
    return I

# Interpolation part

def divided_diff(i, x, f):
    if len(i) == 1:
        return f[0]
    if x[-1] == x[0]:
        return 0
    
    return (divided_diff(i[1:], x[1:], f[1:]) - divided_diff(i[:-1], x[:-1], f[:-1])) / (x[-1] - x[0])


# Tridiagonal matrix algorithm
# See: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
def tridiagonal(a, b, c, d):
    a_res = np.array(a, dtype = np.float64)
    b_res = np.array(b, dtype = np.float64)
    c_res = np.array(c, dtype = np.float64)
    d_res = np.array(d, dtype = np.float64)
    
    for i in range(1, len(d_res)):
        w = a_res[i - 1] / b_res[i - 1]
        b_res[i] = b_res[i] - w * c_res[i - 1]
        d_res[i] = d_res[i] - w * d_res[i - 1]

    x = np.empty(len(d_res), dtype = np.float64)
    x[-1] = d_res[-1] / b_res[-1]
    for i in range(len(d_res) - 2, -1, -1):
        x[i] = (d_res[i] - c_res[i] * x[i + 1]) / b_res[i]

    return x

def spline_interp(nodes, values, x_arr, lambda_0 = 1, mu_n = 1):
    
    diag = np.full(len(nodes), 2)
    mu = np.empty(len(nodes) - 1)
    mu.fill(0)
    for i in range(0, len(nodes) - 2):
        mu[i] = (nodes[i+1] - nodes[i])/(nodes[i+2] - nodes[i])
    mu[-1] = mu_n # Boundary conditions
    
    lambdas = np.empty(len(nodes) - 1)
    lambdas.fill(0)
    lambdas[1:] = 1 - mu[:-1]
    lambdas[0] = lambda_0 # Boundary conditions

    d = np.empty(len(nodes))
    for i in range(0, len(nodes)):
        jj = (max(0, i - 1), i, min(len(nodes) - 1, i + 1))
        d[i] = 6 * divided_diff(jj, np.take(nodes, jj), np.take(values, jj))

    d[0] = d[-1] = 0
            
    M = tridiagonal(mu, diag, lambdas, d)

    res = np.empty(len(x_arr))
    for j, x in enumerate(x_arr):
        i = np.argmax(nodes > x)
        hi = (nodes[i] - nodes[i-1])
        res[j] = M[i-1] * (nodes[i] - x)**3 / (6 * hi)                    \
               + M[i] * (x - nodes[i-1])**3 / (6 * hi)                    \
               + (values[i-1] - M[i-1] * hi**2 / 6) * (nodes[i] - x) / hi \
               + (values[i] - M[i] * hi**2 / 6) * (x - nodes[i-1]) / hi
        
    return res

plt.scatter(x, f_val)
plt.show()

print(f"Trapezoidal             : {trapezoidal(x, f_val)}")
print(f"Simpson                 : {simpson(x, f_val)}")
print(f"Richardson (Trapezoidal): {richardson(x, f_val, trapezoidal, 2)}")
print(f"Richardson (Simpson)    : {richardson(x, f_val, simpson, 4)}")
print(f"Reference (Wolfram)     : {0.367879441171442}")

x_i = np.array([ 0.1,     0.5,      0.9,     1.3,     1.7    ], dtype=np.float64)
f_i = np.array([-2.3026, -0.69315, -0.10536, 0.26236, 0.53063], dtype=np.float64)
k = 80

t = np.linspace(min(x_i), max(x_i), 1001)
plt.plot(t, spline_interp(x_i, f_i, t))
plt.scatter(x_i, f_i)
plt.show()

x = t
f_val = np.sin(k * x) * spline_interp(x_i, f_i, x)

plt.plot(x, f_val)
plt.show()

print(f"Trapezoidal             : {trapezoidal(x, f_val)}")
print(f"Simpson                 : {simpson(x, f_val)}")
print(f"Richardson (Trapezoidal): {richardson(x, f_val, trapezoidal, 2)}")
print(f"Richardson (Simpson)    : {richardson(x, f_val, simpson, 4)}")

