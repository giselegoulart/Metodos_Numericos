import numpy as np

def LagrangeInterp(data, x):
    #Number of data points
    n=len(data)
    #Number of x points
    nx = len(x)

    #Parse x, y data points
    dx = [d[0] for d in data]
    dy = [d[1] for d in data]

    #Allocate space for L(x)
    L = [0.0]*(nx)

    def b(j,xi):
        """Calculate b_j(x_xi)"""
        v = 1.0
        for k in xrange(n):
            if k != j:
                v *= (xi-dx[k]) / (dx[j]-dx[k])
        return v

    #Construct L(x)
    for i,xi in enumerate(x):
        #Construct each element of L(x)
        for j in xrange(n):
            L[i] += dy[j]*b(j,xi)
            
    return L

if '__main__' in __name__:
    import matplotlib as mpl
    mpl.use("TKAgg")
    import matplotlib.pylab as plt
    plt.ion()
    x = lambda n: np.linspace(-1,1,n)
    f = lambda x: np.cos(np.sin(np.pi*x))
    #plt.plot(x(300),f(x(300)),'k')
    n=18
    #LX=x(90)
    x = [1, 2, 5, 6, 7, 8, 10, 13, 17, 20, 23, 24, 25, 27, 27.7, 28, 29, 30]
    y = [3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1, 4.3, 4.1, 3.0]
    
    
    data=zip(x,y)
    LY = LagrangeInterp(data, np.arange(1,30,0.5))
    plt.plot(np.arange(1,30,0.5),LY,'r')