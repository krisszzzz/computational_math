import numpy as np
import matplotlib.pyplot as plt
from math import *


def parse_file(file_path):
    year = []
    population = []
    with open(file_path, 'r') as f:
        lines = f.readlines()

        for line in lines:
            parts = line.strip().split()
            year.append(float(parts[0]))
            population.append(float(parts[1]))

    return year, population 

file_path = 'file.txt'  # Replace 'file.txt' with the path to your file
year, population = parse_file(file_path)

epsilon = 1e-6

target_year = year[5]
target_population = population[5]

year = np.delete(year, (4, 5))
population = np.delete(population, (4, 5))

def show_res(interp_func, x, y):
    t = np.linspace(min(year), max(year), 200)

    plt.plot(t, interp_func(year, population, t))
    plt.plot(target_year,target_population,'ro') 
    plt.scatter(year, population)
    plt.show()

    it = float(interp_func(year, population, x))
    print('Interpolation: {:.0f}'.format(it))
    print('Function: {:.0f}'.format(y))
    print('Error: {:.2f} %'.format(100 * abs(y - it) / y))

def norm(x):
    return np.sqrt(np.dot(x, x))

def divided_diff(i, x, f):
    if len(i) == 1:
        return f[0]
    if x[-1] == x[0]:
        return 0
    
    return (divided_diff(i[1:], x[1:], f[1:]) - divided_diff(i[:-1], x[:-1], f[:-1])) / (x[-1] - x[0])

def newton_interp(nodes, values, x):
    coeffs = np.array([divided_diff(range(i), nodes[:i], values[:i]) for i in range(1, len(nodes) + 1)])
    res = 0
    for i in range(len(coeffs) - 1, 0, -1):
        res = (res + coeffs[i]) * (x - nodes[i - 1])
    return res + coeffs[0]

show_res(newton_interp, target_year, target_population)

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
            
show_res(spline_interp, [target_year], target_population)

def lsqm_interp(nodes, values, x):
    matrix = np.empty((len(nodes), 4))
    for i, node in enumerate(nodes):
        matrix[i] = np.array([1, node, node**2, node**3], dtype = np.float64)
        
    d = np.matmul(matrix.T, values)
    phi = np.matmul(matrix.T, matrix)
    
    coeffs = np.linalg.solve(phi, d)
    
    return coeffs[0] + coeffs[1] * x + coeffs[2] * x**2 + coeffs[3] * x**3


show_res(lsqm_interp, target_year, target_population)
