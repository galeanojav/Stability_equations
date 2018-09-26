#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 10:36:42 2018

@author: javiergaleano
"""

# Sommers's stability

import numpy as np
import matplotlib.pyplot as plt




def random_matrix(S,center,sigma, ty=None):
    """
    Function to calculate the random matrix. 
    S = size of the matrix SxS
    center = value of the matrix diagonal
    sigma = var
    ty [optinal]= None --> normal distibution (0,sigma)
    ty = 'halfnorm' --> half normal distribution
    ty = "exponential" --> exponential distribution scale sigma
    """
    
    ### matrix with "center" in the diagonal
    diago = np.full_like(np.arange(S), center, dtype=np.float32)
    diagonal = np.diagflat(diago)
    
    ### Array with normal distribution mean = 0.0 and var = sigma
    if ty == 'exponential':
        matriz = np.random.exponential(sigma, size=(S*S-S))
        
    if ty == 'halfnorm':
        matriz = np.abs(np.random.normal(loc = 0.0, scale =sigma, size=(S*S-S)))
        
    else:
        matriz = np.random.normal(loc = 0.0, scale =sigma, size=(S*S-S))
    

    ### Random matrix complete
    
    k=0
    for i in range(S):
        for j in range(S):
            if i != j:
                diagonal[i][j]= matriz[k]
                k +=1
    return diagonal

def plot_eig(autovalor):
    """
    Plot eigenvalues. Axis X real part, Axis Y imaginary part
    """
    fig, ax1 = plt.subplots(1, 1)
    
    X1=[]
    Y1=[]
    for i in autovalor:
        X1.append(i.real)
        Y1.append(i.imag)
    ax1.scatter(X1,Y1,alpha=0.3)
    ax1.set(xlabel="Real", ylabel="Im")
    return


S = 1000  # number of species
center = -10.0
sigma = 1/S
diagonal = random_matrix(S,center,sigma, ty='halfnorm')

### mean and var of the matrix

media1 = np.mean(diagonal, dtype = np.float64)
varianza1 = np.var(diagonal,dtype=np.float64)
print(media1, varianza1)


### calculating eigenvalues and plotting

autovalor1 = np.linalg.eigvals(diagonal)

plot_eig(autovalor1)
plt.show()


     
        

#mean, var, skew, kurt = halfnorm.stats(moments='mvsk')

#x = np.linspace(halfnorm.ppf(0.01), halfnorm.ppf(0.99), 100)

#ax.plot(x, halfnorm.pdf(x), 'ro', lw=5, alpha=0.6, label='halfnorm pdf')

#r = halfnorm.rvs(size=10000)
#ax.hist(r, normed=True, histtype='stepfilled', bins = 20, alpha=0.2)
#ax.legend(loc='best', frameon=False)
#plt.show()



