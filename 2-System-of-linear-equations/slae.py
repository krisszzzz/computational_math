
import matplotlib.pyplot as plt;
import sympy as sp;
import numpy as np;


def slae_elem(i, j):
    if ( i == j ):
        return 1
    return 1 / ((i + 1) + (j + 1))

def f_elem(i):
    return 1 / (i + 1)

n = 10

slae_mat = np.fromfunction(np.vectorize(slae_elem, otypes=[float]), (n, n))
f_vec = np.fromfunction(f_elem, (n,))

# solve equation using gauss elimination
def gauss(matrix, f):
    assert len(matrix) == len(f)

    # augmented matrix [A|f]
    augmented_mat = np.hstack((matrix, f.reshape(-1, 1)))   

    # gaussian elimination with partial pivoting 
    for i in range(len(matrix)):

        # find maximum of column ( pivot ) 
        pivot_indx = np.argmax(np.abs(matrix[:,i]))
        # swap pivot to i row 
        augmented_mat[[i, pivot_indx]] = augmented_mat[[pivot_indx, i]] 
        # normalize row with pivot
        pivot = augmented_mat[i][i]
        pivot_row = augmented_mat[i]
        pivot_row /= pivot

        # eliminate the column below the pivot ( excluding the pivot itself )
        for row in augmented_mat[i + 1:]:
            row -= row[i] * pivot_row

    # backward substitution 
    for i in range(len(matrix) - 1, 0, -1):
        row = augmented_mat[i]
        for up in reversed(augmented_mat[:i]):
            up -= up[i] * row

    return augmented_mat[:, -1]

def lu(matrix):
    # assign U = A, L = E
    U = matrix.copy()
    L = np.eye(len(matrix)) 
    # Lji = Uji / Uii 
    for i in range(len(matrix)):
        L[i+1:,i] = U[i+1:,i] / U[i,i]
        U[i+1:] -= L[i+1:,i].reshape(-1, 1) * U[i] 
    return L, U

def lu_solve(matrix, f):
    assert len(matrix) == len(f)

    L, U = lu(matrix)
    x = np.zeros(len(matrix))
    y = np.zeros(len(matrix))
    
    # forward substitution for system Lx = f
    for i in range(len(matrix)):
        x[i] = f[i] - np.dot(L[i], x)

    # backward substitution for system Uy = x
    for i in range(len(matrix) - 1, -1, -1):
        y[i] = (x[i] - np.dot(U[i], y)) / U[i][i]

    return y

def norm(x):
    return np.sqrt(np.dot(x, x))
    
def count_discrepancy(method, matrix, f):
    x = method(matrix, f)
    return norm(np.matmul(matrix, x.transpose()) - f) 


def count_sim_discrepancy(method, max_n, matrix, f):
    discrepancy = np.zeros(max_n)
    for i in range(max_n):
        x = method(matrix, f, i + 1)
        discrepancy[i] = norm(np.matmul(matrix, x.transpose()) - f)
    return discrepancy

# Simple iteration method searching for solution of linear equations
# using formula Xn+1 = BXn + F
def sim(B, F, max_n):
    x = np.zeros(len(F))
    for i in range(max_n):
        x = np.matmul(B, x) + F
    return x
    

def ldu(matrix):
    L = matrix - np.triu(matrix)
    D = np.triu(matrix) + np.tril(matrix) - matrix
    U = matrix - np.tril(matrix)
    # Здесь можно на точное равенство проверить
    return  L, D, U

# jacobi's method
def jacobi_sim(matrix, f, max_n):
    L, D, U = ldu(slae_mat)

    B = -np.matmul(np.linalg.inv(D), (L + U))
    F =  np.matmul(np.linalg.inv(D), f)

    return sim(B, F, max_n)

# gauss-seidel's method
def seidel_sim(matrix, f, max_n):
    L, D, U = ldu(matrix)

    B = -np.matmul(np.linalg.inv(L + D), U)
    F =  np.matmul(np.linalg.inv(L + D), f)

    return sim(B, F, max_n)

# succesive over-relaxion
def sor_sim(matrix, f, max_n, w):
    L, D, U = ldu(matrix)

    B = -np.matmul(np.linalg.inv(D + w * L), (w - 1) * D + w * L)
    F =  np.matmul(np.linalg.inv(D + w * L), f) * w

    return sim(B, F, max_n)

sim_methods = [ jacobi_sim,
                seidel_sim,
                lambda matrix, f, n: sor_sim( matrix, f, n, 1 ),
                lambda matrix, f, n: sor_sim( matrix, f, n, 2 ), 
                lambda matrix, f, n: sor_sim( matrix, f, n, 3 ) ]

print('Discrepancy using gauss elimination:', count_discrepancy(gauss, slae_mat, f_vec))
print('Discrepancy using LU-decomposition:', count_discrepancy(lu_solve, slae_mat, f_vec))

for i in range(len(sim_methods)):
    sim_discrepancy = count_sim_discrepancy(sim_methods[i], 25, slae_mat, f_vec) 
    plt.scatter(range(len(sim_discrepancy)), sim_discrepancy)
    plt.xlabel(r'N')
    plt.ylabel(r'$|\Delta|$')
    plt.show()
