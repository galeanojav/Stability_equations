#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 10:36:42 2018

@author: javiergaleano
"""

# Sommers's stability

import numpy as np
import matplotlib.pyplot as plt


S = 1000  # number of species

### Sommers et al. Stability

vari = 1/S

somer = np.zeros((S,S))
ma_som = np.random.normal(loc = 0.0, scale =vari, size=(S*S-S))

t=0
for i in range(S):
    for j in range(S):
        if i != j:
            somer[i][j]= ma_som[t]
            t +=1
            
### mean and var of the matrix

mediasom = np.mean(somer, dtype = np.float64)
varianzasom = np.var(somer,dtype=np.float64)

print(mediasom, varianzasom)

### calculating eigenvalues and plotting

autovalor = np.linalg.eigvals(somer)

X=[]
Y=[]
for i in autovalor:    
    X.append(i.real)
    Y.append(i.imag)


fig, ax = plt.subplots(1, 1)
ax.scatter(X,Y,alpha=0.3)
ax.grid(True)
ax.set_title("Sommers'stability")
ax.set(xlabel="Real", ylabel="Im")



#mean, var, skew, kurt = halfnorm.stats(moments='mvsk')

#x = np.linspace(halfnorm.ppf(0.01), halfnorm.ppf(0.99), 100)

#ax.plot(x, halfnorm.pdf(x), 'ro', lw=5, alpha=0.6, label='halfnorm pdf')

#r = halfnorm.rvs(size=10000)
#ax.hist(r, normed=True, histtype='stepfilled', bins = 20, alpha=0.2)
#ax.legend(loc='best', frameon=False)
#plt.show()



