# -*- coding: utf-8 -*-

# Exercicio 3

import numpy as np
import matplotlib.pylab as pl
#==============================================================================

x = [0, 1, 2, 6, 7]
y = [0, 0, 1, 2, 3]
n = len(x)

comb = []
for i in range(n):
    for j in range(i+1, n):
       comb.append([i,j])

a = []  
b = []
# combinacao de 2 em 2 pontos
for p in comb:
    m = (y[p[1]] - y[p[0]])/float(x[p[1]] - x[p[0]])
    a.append(m)
    t = -m*x[p[0]]+y[p[0]]
    b.append(t)

a_par = []
b_par = []
dist = []
for i in range(len(a)):
    for j in range(i+1, len(a)):
        if(a[i]==a[j] and a[i]!=0):
           a_par.append(a[i])
           b_par.append([b[i], b[j]])

# calculo do erro (distancias) entre os pontos e as retas com mesmo coef angular
# Distancia maxima selecionada como erro => norma infinita
for l in range(len(a_par)):
    distance = 0
    for i in range(n):
        d1 =  abs(a_par[l]*x[i] + y[i] + b_par[l][0])/float(np.sqrt((a_par[l])**2 + 1))
        d2 =  abs(a_par[l]*x[i] + y[i] + b_par[l][1])/float(np.sqrt((a_par[l])**2 + 1))
        print a_par[l], d1, d2
        d = max([d1,d2])
        # Seleciona maxima distancia
        if(d>distance):
            distance = d
    dist.append(distance)

# menor distancia entre as maximas  
min1 = min(dist)
min1_pos = dist.index(min1)
a1 = a_par[min1_pos]
b1 = b_par[min1_pos]

# media do b das duas retas com melhor coeficiente angular
b_mean = (b1[0]+b1[1])/2

# Plot do grafico
y_adjust1 = [(a1*xi+b1[0]) for xi in x]
y_adjust2 = [(a1*xi+b1[1]) for xi in x]
y_adjust = [(a1*xi+b_mean) for xi in x]
pl.plot(x,y_adjust1)
pl.plot(x,y_adjust2)
pl.plot(x,y_adjust)
pl.title('Erro: '+str('%.2f'% min1)+' - Coef. Angular: '+str(a1))
pl.xlabel('x')
pl.ylabel('y')
pl.scatter(x,y)

