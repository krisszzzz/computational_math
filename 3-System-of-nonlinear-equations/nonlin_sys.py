
import matplotlib.pyplot as plt;
import numpy as np;
import sympy as sp;

sol = [ 0.51015, -0.201838 ]

def norm(x):
    return np.sqrt(np.dot(x, x))
 
def residual_plot(f_iter, iter_count = 10):
    x = np.zeros((iter_count, 2))
    x[0] = [0.9, 0] # points that are in localization area
    for i in range(len(x) - 1):
        x[i + 1] = f_iter(x[i])

    plt.scatter(range(len(x)), [ norm(x[i] - sol) for i in range(len(x)) ])
    plt.xlabel('iteration number')
    plt.ylabel('$\Delta x$')
    plt.show()

    print('solution [x, y] = {}',format(x[-1]))

# Iterative method
f_iter = lambda x: [ 1 - 0.5 * np.cos(x[1]),
                         np.sin( x[0] + 1 ) - 1.2 ]

residual_plot(f_iter)

# Newton method
f = lambda x: [ np.sin(x[0] + 1) - x[1] - 1.2,
                    2*x[0] + np.cos(x[1]) - 2 ]

J_inv = lambda x: 1/(np.cos(x[0] + 1) * np.sin(x[1]) -2) * np.array([[np.sin(x[1]),      -1               ],
                                                                     [2,                -np.cos(x[0] + 1)]])
f_newton = lambda x: x - np.matmul(J_inv(x), f(x))

residual_plot(f_newton)
