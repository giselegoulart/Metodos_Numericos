# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

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

        # LaTeX: u_{ij} = a_{ij} - \sum_{k=1}^{i-1} u_{kj} l_{ik}                                                                                                                                                                                      
        for i in range(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in xrange(i))
            U[i][j] = PA[i][j] - s1

        # LaTeX: l_{ij} = \frac{1}{u_{jj}} (a_{ij} - \sum_{k=1}^{j-1} u_{kj} l_{ik} )                                                                                                                                                                  
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
    
def read_file(name_input):
    try:
      fid=open('./'+name_input,'r')
      aux = pd.read_table(fid, sep='\t', header=None)
      fid.close()
    except EOFError:
      print 'Erro na leitura do arquivo'
    except IOError:
      print 'Arquivo inexistente'
    
    #col_names = ['tempo', '\n']
    #aux.columns = col_names
    return pd.DataFrame(aux)

#==============================================================================

files_A = ['A-4.dat', 'A-6.dat', 'A-8.dat', 'A-10.dat']
files_b = ['b-4.dat', 'b-6.dat', 'b-8.dat', 'b-10.dat']

for arq in range(len(files_A)):
    A = read_file(files_A[arq])
    del A[len(A)]
    b = read_file(files_b[arq])
    A = np.array(A)
    b = np.array(b)
    n = len(b)
    print 20*'#', 'Sistema ', n, 'x', n, 20*'#'
    #print A, b

    P, L, U = lu_decomposition(A)
    
    # Solucao com substituicao
    ysolution = lowerSubs(L,b)
    
    xsolution = upperSubs(U,ysolution)
    print 'Solucao com rotina'
    print xsolution
    
    # Solucao com Numpy
    solution_np = np.linalg.solve(L,b)
    solution_np = np.linalg.solve(U,solution_np)
    solution_np = solution_np.reshape(1,n)
    print 'Solucao Numpy com A=LU'
    print solution_np
    
    print 'Solucao numpy com A'
    solution_np_A = np.linalg.solve(A,b)
    solution_np_A = solution_np_A.reshape(1,n)
    print(solution_np_A), '\n'
    
    # Calculo da matriz inversa
    AI = []
    exact_solution = []
    for col in range(n):
        ej=[]
        for i in range(n):
            if(i==col):
                ej.append(1)
            else:
                ej.append(0)
        bj = lowerSubs(L,ej)  
        bj = upperSubs(U,bj)
        AI.append(bj)
        exact_solution.append(101)
        
    AI = np.array(AI)
    AI = AI.transpose()

    AI_norm = np.linalg.norm(AI, np.inf)
    A_norm = np.linalg.norm(A, np.inf)
    cond_A = AI_norm*A_norm
    print '\nCondicionamento de A: ', cond_A, '\n'
       
    xsolution_norm = np.linalg.norm(xsolution, np.inf)  
    x_exact_norm = np.linalg.norm(exact_solution, np.inf)
    dif_norm = abs(xsolution_norm-x_exact_norm)
    print 'Diferenca em norma infinito entre a solucao esperada e a solucao aproximada: '
    print dif_norm, '\n'
