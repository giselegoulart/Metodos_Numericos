from math import *
import matplotlib.pylab as pl
import numpy as np

#==============================================================================

def lagrange_interpolation(x, y, u):
  r = range(len(y))
  a =       [y[i]/product( x[i]-x[j] for j in r if j != i ) for i in r]
  return sum(a[i]*product([u   -x[j] for j in r if j != i]) for i in r)


def product(a): 
  p = 1
  for i in a: p *= i
  return p
  
#==============================================================================
# Entrada de dados  
x = [1, 2, 5, 6, 7, 8, 10, 13, 17, 20, 23, 24, 25, 27, 27.7, 28, 29, 30]
y = [3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1, 4.3, 4.1, 3.0]


"""
Polinomio de grau n-1 (Lagrange)
"""
approx_total = []
for i in np.arange(1,30,0.2):
    estim = lagrange_interpolation(x, y, i)
    approx_total.append(estim)


pl.ylim(-70,110)
pl.grid(True)
pl.title("Aproximacao de polinomio de grau "+str(len(x)-1))
pl.xlabel('x')
pl.ylabel('y')
pl.plot(np.arange(1,30,0.2), approx_total) 

"""
Polinomios por partes (Lagrange)
"""

# Divisao dos intervalos para aproximacao por partes das funcoes
# Sao utilizadas funcoes de primeiro e segundo grau

aux = []
for i in np.arange(1,2, 0.2):
    x_lin = [x[0], x[1]]
    y_lin = [y[0], y[1]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim)
    
for i in np.arange(3,5, 0.2):
    x_lin = [x[1], x[2]]
    y_lin = [y[1], y[2]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim)
    
for i in np.arange(6,12, 0.2):
    x_lin = [x[3],x[5],x[6]]
    y_lin = [y[3],y[5],y[6]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim)
    
for i in np.arange(13.4,17, 0.2):
    x_lin = [x[6],x[8]]
    y_lin = [y[6],y[8]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim)
    
for i in np.arange(18,23, 0.2):
    x_lin = [x[8],x[9],x[10]]
    y_lin = [y[8],y[9],y[10]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim)
    
for i in np.arange(24,25, 0.2):
    x_lin = [x[11],x[12]]
    y_lin = [y[11],y[12]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim)
    
for i in np.arange(25.5,28, 0.2):
    x_lin = [x[12],x[13],x[15]]
    y_lin = [y[12],y[13],y[15]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim)
    
for i in np.arange(28.5,30, 0.2):
    x_lin = [x[15],x[16],x[17]]
    y_lin = [y[15],y[16],y[17]]
    estim = lagrange_interpolation(x_lin, y_lin, i)
    aux.append(estim) 
    
pl.figure()
pl.plot(np.linspace(1,30,len(aux)), aux,'b-')
#pl.plot(x,y,'ro')
pl.title("Aproximacao por partes com polinomios de grau 1 e 2")
pl.xlabel('x')
pl.ylabel('y')
pl.ylim(0,20 )
plt.grid(True)

#pl.close('all')

"""
Ajuste de Spline Natural
"""

n = len(x)
h = [x[i+1]-x[i] for i in range(n-1)]


b = [(y[i+1]-y[i])/float(h[i]) for i in range(n-1)]

m = len(h)

# Para determinar as incognitas z um sistema de equacoes eh montado
# Para splines naturais z0 = zn = 0
# Eliminando a primeira e ultima coluna temos matriz quadrada (n-1)x(n-1)

s = (n-2,n-1)
p = n-2
q = n-1
A = np.zeros(s)

if s == (1,2):
    A[0][0]=2*(h[0]+h[1])
    A[0][1]=6*(b[1]-b[0])
    
if s == (2,3):
    A[0][0]=2*(h[0]+h[1])
    A[0][1]=h[1]
    A[1][0]=h[1]
    A[1][1]=2*(h[1]+h[2])
    
    A[0][2]=6*(b[1]-b[0])
    A[1][2]=6*(b[2]-b[1])

if (n-2) > 2:
    A[0][0]=2*(h[0]+h[1])
    A[0][1]=h[1]
    A[p-1][q-3]=h[m-2]
    A[p-1][q-2]=2*(h[m-2]+h[m-1])
    
    A[0][q-1]=6*(b[1]-b[0])
    A[p-1][q-1]=6*(b[m-2]-b[m-3])
    
    for i in range(1,p-1):
        A[i][i-1]= h[i]
        A[i][i]= 2*(h[i]+h[i+1])
        A[i][i+1]=h[i+1]
        
        A[i][q-1] = 6*(b[i+1]-b[i])
        

# Matriz aumentada
# Solucao do sistema por eliminacao gaussiana
(r,t) = np.shape(A)

if(r,t) == (1,2):
    z = A[0][1]/A[0][0]

if r > 1:
    for i in range(0,r):
        pivot = A[i][i]
    
        for j in range(i+1,r):
            l = (A[j][i]/pivot)
        
            for k in range(i,t):
                A[j][k] = A[j][k]-(l*A[i][k])

    z = A[:,t-1]

    for i in range(r-1,-1,-1):
        for j in range(r-1,i,-1):
            z[i] = z[i]-(A[i][j]*z[j])

        z[i]=z[i]/A[i][i]
        A[i][t-1] = z[i]

z = np.append(np.array([0]),z)
z = np.append(z,0)


# Calculo dos coeficientes a[i], b[i], c[i] e d[i]
n = len(x)-1

e = np.zeros(n)
f = np.zeros(n)
g = np.zeros(n)
o = np.zeros(n)

for i in range(0,n-1):
    e[i]=(z[i+1]-z[i])/(6*h[i])
    f[i]=(z[i]/2)
    g[i]=b[i]-h[i]*((z[i+1]+2*z[i])/6)
    o[i]=y[i]

o[n-1]=y[n-1]
g[n-1]=3*e[n-2]*(h[n-2]**2)+2*f[n-2]*(h[n-2])+g[n-2]
f[n-1]=z[n-1]/2
e[n-1]=((y[n]-y[n-1])-(f[n-1]*(h[n-1]**2)+g[n-1]*h[n-1]))/(h[n-1]**3)


# Plot da spline polinomial
n = len(x)-1

t = []
u = []
axis = []

y_aprox = []
x_disc = []
from scipy.integrate import simps, romb
for i in range(0,n):
    t = np.linspace(x[i],x[i+1],100)
    axis = np.append(axis,t)
    u = np.append(u,e[i]*((t-x[i])**3)+f[i]*((t-x[i])**2)+g[i]*(t-x[i])+o[i])
    y_aprox.append(e[i]*((t-x[i])**3)+f[i]*((t-x[i])**2)+g[i]*(t-x[i])+o[i])
    x_disc.append(t)
    
pl.figure()
pl.plot(axis, u, 'b-')
#pl.plot(x,y,'ro')
pl.title("Aproximacao por spline natural")
pl.xlabel('x')
pl.ylabel('y')
pl.grid(True)
pl.ylim(0,20)
pl.show()


"""
Integracao da area abaixo da curva determinada pelas splines com 1/3 de Simpson
"""

def integrate(y_vals, h):
    i = 1
    total = y_vals[0] + y_vals[-1]
    for y in y_vals[1:-1]:
        if i % 2 == 0:
            total += 2 * y
        else:
            total += 4 * y
        i += 1
    return total * (h / 3.0)
    


area=0
area1=0
for i in range(len(x_disc)):
    h = float(x_disc[i][1]-x_disc[i][0])
    area1 += integrate(y_aprox[i], h)
    #funcao do scipy
    #area+=simps(y_aprox[i], x_disc[i])
print "Area abaixo da curva: ", area1