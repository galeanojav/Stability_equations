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
sigma = vari
### matrix with -1.0 in the diagonal
diago = np.full_like(np.arange(S), -10.0, dtype=np.float32)
diagonal1 = np.diagflat(diago)
diagonal2 = np.diagflat(diago)
diagonal3 = np.diagflat(diago)
diagonal4 = np.diagflat(diago)


### Array with normal distribution mean = 0.0 and var = sigma
    
matriz1 = np.random.normal(loc = 0.0, scale =sigma, size=(S*S-S))
matriz2 = np.random.triangular(-20, -10, 3, size=(S*S-S))
matriz3 = np.random.wald(0.2, sigma, size=(S*S-S))
matriz4 = np.random.exponential(sigma, size=(S*S-S))

### Random matrix complete
    
k=0
for i in range(S):
    for j in range(S):
        if i != j:
            diagonal1[i][j]= matriz1[k]
            diagonal2[i][j]= matriz2[k]
            diagonal3[i][j]= matriz3[k]
            diagonal4[i][j]= matriz4[k]
            k +=1
       

### mean and var of the matrix

media1 = np.mean(diagonal1, dtype = np.float64)
varianza1 = np.var(diagonal1,dtype=np.float64)

print(media1, varianza1)

media2 = np.mean(diagonal2, dtype = np.float64)
varianza2 = np.var(diagonal2,dtype=np.float64)

print(media2, varianza2)

media3 = np.mean(diagonal3, dtype = np.float64)
varianza3 = np.var(diagonal3,dtype=np.float64)

print(media3, varianza3)

media4 = np.mean(diagonal4, dtype = np.float64)
varianza4 = np.var(diagonal4,dtype=np.float64)

print(media4, varianza4)

### calculating eigenvalues and plotting

autovalor1 = np.linalg.eigvals(diagonal1)
autovalor2 = np.linalg.eigvals(diagonal2)
autovalor3 = np.linalg.eigvals(diagonal3)
autovalor4 = np.linalg.eigvals(diagonal4)

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2)
    
X1=[]
Y1=[]
for i in autovalor1:    
    X1.append(i.real)
    Y1.append(i.imag)

X2=[]
Y2=[]
for j in autovalor2:    
    X2.append(j.real)
    Y2.append(j.imag)
    
X3=[]
Y3=[]
for k in autovalor3:    
    X3.append(k.real)
    Y3.append(k.imag)
    
X4=[]
Y4=[]
for l in autovalor4:    
    X4.append(l.real)
    Y4.append(l.imag)
    
ax1.scatter(X1,Y1,alpha=0.3)
ax1.set(xlabel="Real", ylabel="Im")
    
ax2.scatter(X2,Y2,alpha=0.3)
ax2.set(xlabel="Real", ylabel="Im")
    
ax3.scatter(X3,Y3,alpha=0.3)
ax3.set(xlabel="Real", ylabel="Im")
    
ax4.scatter(X4,Y4,alpha=0.3)
ax4.set(xlabel="Real", ylabel="Im")
    




     
        

#mean, var, skew, kurt = halfnorm.stats(moments='mvsk')

#x = np.linspace(halfnorm.ppf(0.01), halfnorm.ppf(0.99), 100)

#ax.plot(x, halfnorm.pdf(x), 'ro', lw=5, alpha=0.6, label='halfnorm pdf')

#r = halfnorm.rvs(size=10000)
#ax.hist(r, normed=True, histtype='stepfilled', bins = 20, alpha=0.2)
#ax.legend(loc='best', frameon=False)
#plt.show()



