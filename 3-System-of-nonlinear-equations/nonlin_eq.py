
import matplotlib.pyplot as plt;
import numpy as np;
import sympy as sp;

x = sp.Symbol('x')
f = x*2**x - 1

sol = 0.6411857445 # approximate solution of equation

def residual_plot(f_iter, iter_count = 10):
    x = np.zeros(iter_count) # array of X_k
    x[0] = 0.6 # start position that belongs to [0,5; 1] 

    for i in range(len(x) - 1):
        x[i + 1] = f_iter(x[i])

    plt.scatter(range(len(x)), abs(x - sol))
    plt.xlabel('iteration number')
    plt.ylabel('$\Delta x$')
    plt.show()

    print('solution [x] = {}',format(x[-1]))

# Iterative method
f_iter = lambda x: 1 / (2**x)
residual_plot(f_iter)
# Newton method
f_newton = sp.lambdify( x, x - f / f.diff() ) 
residual_plot(f_newton, 5)
