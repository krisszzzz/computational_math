
import matplotlib.pyplot as plt;
import sympy as sp;
import numpy as np;

x = sp.Symbol('x')

# symbol description of each functions
functions_list = [
    sp.sin(x**2),
    sp.cos(sp.sin(x)),
    sp.exp(sp.sin(sp.cos(x))),
    sp.log(x + 3),
    (x + 3)**0.5
]

# callable derivative of the functions
derivatives = [ sp.lambdify( x, fn.diff()) for fn in functions_list ] 

# callable functions
functions = [ sp.lambdify(x, fn) for fn in functions_list ]

numeric_derivatives = [
    lambda f, x, h: (f(x + h) - f(x)) / h,
    lambda f, x, h: (f(x) - f(x - h)) / h,
    lambda f, x, h: (f(x + h) - f(x - h)) / (2 * h),
    lambda f, x, h: (4 / 3) * (f(x + h) - f(x - h)) / (2 * h) 
                    - (1 / 3) * (f(x + 2 * h) - f(x - 2 * h)) / (4 * h),
    lambda f, x, h: (3 / 2) * (f(x + h) - f(x - h)) / (2 * h)
                    - (3 / 5) * (f(x + 2 * h) - f(x - 2 * h)) / (4 * h)
                    + (1 / 10) * (f(x + 3 * h) - f(x - 3 * h)) / (6 * h),
]

markers = [ 'o', 'v', 's', 'h', 'x' ]

x0 = 1.5
n_max = 21
steps = np.array([(2 / 2**n) for n in range(1, n_max + 1)])

suite = zip(functions_list, derivatives, functions)

for fn_name, deriv, fn in suite:
    for i, num_deriv in enumerate(numeric_derivatives, start = 0):
        diff = np.abs(num_deriv(fn, x0, steps) - deriv(x0))
        plt.plot(np.log(steps), np.log(diff), label = 'Method {}'.format(i + 1), marker = markers[i])  
    
    plt.xlabel(r'$\ln(h)$')
    plt.ylabel(r'$\ln(\Delta)$')
    plt.title('Error for derivative of {} at x = {}'.format(fn_name, x0))
    plt.legend()
    plt.show()
