#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import scipy.linalg as sla
from scipy.optimize import broyden1
 
def broyden_good(x, y, w1, w2, f_equations, J_equations, tol=10e-10, maxIters=50):
    steps_taken = 0
 
    f = f_equations(x,y,w1,w2)
    J = J_equations()
 
    while np.linalg.norm(f,2) > tol and steps_taken < maxIters:
        if np.linalg.det(J)==0:
            steps_taken=maxIters
            break
        s = sla.solve(J,-1*f)
 
        x = x + s[0]
        y = y + s[1]
        newf = f_equations(x,y,w1,w2)
        z = newf - f
 
        J = J + (np.outer ((z - np.dot(J,s)),s)) / (np.dot(s,s))
 
        f = newf
        steps_taken += 1
 
    return steps_taken, x, y, w1, w2  
    
    
'''
Metodo de Broyden com modificacoes para a solucao do trabalho 7 de quadratura de Gauss
'''    
def broyden_good_gauss(x, y, w1, w2, w3, w4, f_equations, J_equations, tol=10e-10, maxIters=50):
    steps_taken = 0
 
    f = f_equations(x,y,w1,w2, w3, w4)
    J = J_equations()
 
    while np.linalg.norm(f,2) > tol and steps_taken < maxIters:
 
        s = sla.solve(J,-1*f)
 
        x = x + s[0]
        y = y + s[1]
        newf = f_equations(x,y,w1,w2, w3, w4)
        z = newf - f
 
        J = J + (np.outer ((z - np.dot(J,s)),s)) / (np.dot(s,s))
 
        f = newf
        steps_taken += 1
 
    return steps_taken, x, y, w1, w2, w3, w4  
 
'''
Sistemas
''' 
def system1_3(xi,yi,zi,wi):
    return np.array([(xi**2+yi**2+zi**2-1)*yi**2, (2*xi**2+yi**2-4*zi)*yi**2, (3*xi**2-4*yi+zi**2)*yi**2])
    
def system2_3(xi,yi,zi,wi):
    return np.array([xi+xi**2-2*yi*zi-0.1, yi-yi**2+3*xi*yi+0.2, zi+zi**2+2*xi*yi-0.3])

def system3_3(xi,yi,zi,wi):
    return np.array([10*xi-2*yi**2+yi-2*zi-5, 8*yi**2+4*zi**2-9, 8*yi*zi+4])

def system4_3(xi,yi,zi,wi):
    return np.array([xi**2+yi-37, xi-yi**2-5, xi+yi+zi-3])
    
def system5_3(xi,yi,zi,wi):
    return np.array([xi**2+2*yi**2-yi-2*zi, xi**2-8*yi**2+10*zi, xi**2/(7*yi*zi)-1])

def system6_3(xi,yi,zi,wi):
    return np.array([3*xi-np.cos(yi*zi)-1/2, xi**2-81*(yi+0.1)**2+np.sin(zi)+1.06, np.e**(-xi*yi)+20*zi+(10*np.pi-3)/3])

def system7_3(xi,yi,zi,wi):
    return np.array([xi+np.cos(xi*yi*zi)-1, (1-xi)**(1/4)+yi+0.05*zi**2-0.15*zi-1, -xi**2-0.1*yi**2+0.01*yi+zi-1])

def system8_4(xi,yi,zi,wi):
    return np.array([xi+10*yi, (yi-zi)**2, np.sqrt(5)*(zi-wi), np.sqrt(10)*(xi-wi)**2])

def system9_4(x1,x2,w1,w2):
    return np.array([w1+w2-2, w1*x1+w2*x2, w1*x1**2+w2*x2**2-2/3, w1*x1**3+w2*x2**3])

def system10_4(x1,x2,x3,x4):
    return np.array([4*x1-x2+x3-x1*x4, -x1+3*x2-2*x3-x2*x4, x1-2*x2+3*x3-x3*x4, x1**2+x2**2+x3**2-1])
    
def system_quad_gauss(t1,t2, w0,w1,w2,w3):
    return np.array([-w0+w1*t1**5+w2*t2**5+w3, w0+w1*t1**4+w2*t2**4+w3-2/5, -w0+w1*t1**3+w2*t2**3+w3, w0+w1*t1**2+w2*t2**2+w3-2/3, -w0+w1*t1+w2*t2+w3, w0+w1+w2+w3-2])

'''
Matrizes Jacobianas Inciais
''' 
def J2():
    return np.array([[1,0],
             [0, 1]])

def J3():
    return np.array([[1,0,0],
             [0, 1,0], [0,0,1]])
'''
Modificacao da matriz Jacobiana inicial pra ser preenchida randomicamente 
'''             
def J3_R():
    return np.array([[uniform(0,1),uniform(0,1),uniform(0,1)],
             [uniform(0,1), uniform(-0,1),uniform(-0,1)], [uniform(0,1),uniform(0,1),uniform(0,1)]])
	
def J4():
    return np.array([[1,0,0,0],
             [0, 1,0,0], [0,0,1,0],[0,0,0,1]])

'''
Modificacao da matriz Jacobiana inicial pra ser preenchida randomicamente 
'''             
def J4_R():
    return np.array([[uniform(0,1),uniform(0,1),uniform(0,1),uniform(0,1)],
             [uniform(0,1), uniform(0,1),uniform(0,1),uniform(0,1)], [uniform(0,1),uniform(0,1),uniform(0,1),uniform(0,1)],[uniform(0,1),uniform(0,1),uniform(0,1),uniform(0,1)]])

def J6():
    return np.array([[1,0,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,0,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]])
 
 
def f8(x):
    return np.array([x[0]+10*x[1], (x[1]-x[2])**2, np.sqrt(5)*(x[2]-x[3]) , np.sqrt(10)*(x[0]-x[3])**2])

def f2(x):
    return np.array([x[0]+x[0]**2-2*x[1]*x[2]-0.1, x[1]-x[1]**2+3*x[0]*x[1]+0.2, x[2]+x[2]**2+2*x[0]*x[1]-0.3])

def f6(x):
    return np.array([3*x[0]-np.cos(x[1]*x[2])-1/2, x[0]**2-81*(x[1]+0.1)**2+np.sin(x[2])+1.06, np.e**(-x[0]*x[1])+20*x[2]+(10*np.pi-3)/3])

def f9(x):
        return np.array([x[0]+x[1]-2, x[0]*x[2]+x[1]*x[3], x[0]*x[2]**2+x[1]*x[3]**2-2/3, x[1]*x[2]**3+x[1]*x[3]**3])
    
def f10(x):
        return np.array([4*x[0]-[1]+x[2]-x[0]*x[3], -x[0]+3*x[1]-2*x[2]-x[1]*x[3], x[0]-2*x[1]+3*x[2]-x[2]*x[3], x[0]**2+x[1]**2+x[2]**2-1])


#==============================================================================

tol = 10.0**(-10)
maxIters =50000
x0 = 0.860876501058
y0 = -0.440634731817
z0 = 0.570445071156
w0 = 0.980561718853
  
#x0 = uniform(-1,1)
#y0 = uniform(-1,1)
#z0 = uniform(-1,1)
#w0 = uniform(-1,1)
#print x0, y0, z0, w0

systems=[system1_3,system2_3, system3_3, system4_3, system5_3, system6_3, system7_3, system8_4, system9_4, system10_4,]

sol = broyden1(f9,[x0,y0,z0,w0])

print sol

n, x, y, z, w = broyden_good(x0, y0, z0,w0,systems[9] , J4, tol, maxIters)
print n, x, y, z, w
if n == maxIters:
    rod = 0
    while n==maxIters and rod<200:
#        x0 = uniform(-1,1)
#        y0 = uniform(-1,1)
#        z0 = uniform(-1,1)
#        w0 = uniform(-1,1)
        n, x, y, z, w = broyden_good(x0, y0, z0,w0,systems[9] , J4_R, tol, maxIters)
        rod+=1
        print 'Rodada: ', rod
print n, rod, x, y, z, w


'''
Sistema Nao linear do problema de quadratura de gauss
'''
#n,t1,t2,w1,w2,w3,w4 =  broyden_good_gauss(0.1,0.1,0.1,0.1,0.1,0.1,system_quad_gauss , J6, tol, maxIters=500000000)
#print n,t1,t2,w1,w2,w3,w4


