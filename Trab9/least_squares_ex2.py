# -*- coding: utf-8 -*-

# Exercicio 2

import numpy as np
import matplotlib.pylab as pl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#==============================================================================

def mult_matrix(M, N):                                                                                                                                                                                                  
    n = len(M)
    matrizR=[]

    # Multiplicacao das matrizes M, N
    for i in range(n):
        matrizR.append([])
        for j in range(n):
            val = 0
            for k in range(n):
                    val += M[i][k]*N[k][j]
            matrizR[i].append(val)
            
    return matrizR
                                                                                                                                                                                    

def pivot_matrix(M):
    m = len(M)

    # Matriz identidade                                                                                                                                                                                            
    id_mat = [[float(i ==j) for i in xrange(m)] for j in xrange(m)]
    
    # Organizacao da matriz identidade para que o maior elemento de 
    # cada coluna de M fique na diagonal de M                                                                                                                                                                                           
    for j in range(m):
        row = max(range(j, m), key=lambda i: abs(M[i][j]))
        if j != row:
            # Troca de linhas                                                                                                                                                                                                                            
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]

    return id_mat

def lu_decomposition(A):
    n = len(A)

    # Criacao das matrizes L e U zeradas                                                                                                                                                                                                                
    L = [[0.0] * n for i in xrange(n)]
    U = [[0.0] * n for i in xrange(n)]

    # Matriz de pivoteamento                                                                                                                                                                                             
    P = pivot_matrix(A)
    PA = mult_matrix(P, A)

    # Decomposicao LU                                                                                                                                                                                                                    
    for j in range(n):
        # Diagonal de L=1                                                                                                                                                                                                   
        L[j][j] = 1.0
                                                                                                                                                                                      
        for i in range(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in xrange(i))
            U[i][j] = PA[i][j] - s1
                                                                                                                                                                  
        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in xrange(j))
            L[i][j] = (PA[i][j] - s2) / U[j][j]

    return (P, L, U)

def lowerSubs(L, b):
    n = len(L)
    #  Resolve sistema Ly = b
    y = [0 for i in range(n)]
    
    for i in range(0,n,1):
        y[i] = b[i]/float(L[i][i])
        for k in range(0,i,1):
            y[i] -= y[k]*L[i][k]
    return y
 
     
def upperSubs(U, y):
    # Resolve sistema Ux = y
    n = len(U)
    x = np.zeros(n)

    x[-1] = y[-1]/U[-1][ -1] 
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.sum(U[i][i+1:] * x[i+1:]))/U[i][i] 

    return x


def linRegMod(x, f, g):
    n = len(x)
    
    sumx = np.sum(x)
    sumx2 = [i**2 for i in x]
    sumx2 = np.sum(sumx2)
    sumf = np.sum(f)
    sumg = np.sum(g)
    sumfg = [(f[i]+g[i])*x[i] for i in range(n)]
    sumfg = np.sum(sumfg)
    
    A = [[2*sumx2, sumx, sumx],[sumx, n, 0],[sumx, 0, n]]
    b = [sumfg, sumf, sumg]
    
    # Solucao de martriz 3x3 para encontrar coeficientes a, b, c
    P, L, U = lu_decomposition(A)
    
    # Solucao com substituicao
    ysolution = lowerSubs(L,b)
    
    xsolution = upperSubs(U,ysolution)
    
    return xsolution[0], xsolution[1], xsolution[2]

#==============================================================================

x = [0, 1, 2, 3]
f = [-0.2, 0.4, 0.2, 0.6]
g = [1.2, 1.6, 1.8, 2]
n = len(x)

a, b, c = linRegMod(x, f, g)

print a, b, c
y_adjust_f = [a*xi+b for xi in x]
y_adjust_g = [a*xi+c for xi in x]

pl.scatter(x,f)
pl.scatter(x,g)
pl.plot(x, y_adjust_f)
pl.plot(x, y_adjust_g)
pl.legend(['f', 'g'], loc='best')
pl.xlabel('x')
pl.ylabel('y')
pl.show()
