# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import factorial
import matplotlib.pylab as pl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

x_vet =[-20, -10, -2, -1, 1, 2, 10, 20]
erro_vet = []

for x in x_vet:
    soma = 1
    for i in range(1,41):
        soma += x**i/factorial(i, exact=False)
    exact = np.e**x
    erro = abs(soma-exact)
    erro_vet.append(erro)
    print erro

pl.plot(x_vet, erro_vet, '-o')
pl.xlabel('Valores de x')
pl.ylabel('Erro')
pl.title('Erro de Aproximacao de f = e^{x}')
pl.show()
