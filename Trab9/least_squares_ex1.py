# -*- coding: utf-8 -*-

# Exercise 1
# adjust a curve y=ae^(bx)
# linearization lny = lna + bx
# Y = lny
# a0 = lna
# a1 = b
# X = x
# Y = a0 + a1X

import numpy as np
import matplotlib.pylab as pl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#==============================================================================

def linReg(x, y):
    "Least Squares implementation for y=ax + b"
    n = len(x)
    
    sumx = np.sum(x)
    sumx2 = [i**2 for i in x]
    sumx2 = np.sum(sumx2)
    sumxy = [x[i]*y[i] for i in range(n)]
    sumxy = np.sum(sumxy)
    sumy = np.sum(y)
    sumy2 = [i**2 for i in y]
    sumy2 = np.sum(sumy2)
    
    # angular coef a
    a = (sumx*sumy - n*sumxy)/(sumx**2 - n*sumx2)
    b = (sumy - a*sumx)/n
    
    return b, a
    
    
    
#==============================================================================
x = [0, 1, 2, 3, 4, 5 , 6]
y = [32, 47, 65, 92, 132, 190, 275]

y_ln = [np.log(i) for i in y]

a0, a1 = linReg(x, y_ln)

a = np.e**a0
b = a1

y_adjust = [a*np.e**(b*xi) for xi in x]

pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.title("y = a e^{bx}")
pl.xlabel('Horas')
pl.ylabel('Bacterias')
pl.scatter(x, y)
pl.plot(x, y_adjust)
pl.show()

x_2000 = (np.log(2000) - np.log(a))/b
print "Para o numero de bacterias por unidade de volume chegar em 2000: ", "%.2f" % x_2000, " horas sao necessarias."
print "horas >", "%.2f" % x_2000, " numero de bacterias ultrapassa 2000"

    
    
        
