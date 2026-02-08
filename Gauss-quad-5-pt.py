# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 21:38:33 2023

@author: Fire2002
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as spec
import scipy.integrate as sp

'''
List the first five Legendre polynomials,
P_k(x), k = 1,...,5 over [-1,1]
'''

x = np.linspace(-1,1)

for k in range(1,6):
    y = spec.eval_legendre(k,x)
    plt.plot(x,y,label=r'$P_{}(x)$'.format(k))

plt.title("Legendre Polynomials")
plt.xlabel("x")
plt.ylabel(r'$P_k(x)$')
plt.legend(loc='lower right')
plt.show()


def test_function(x):
    # a 72 coeff rules for a complex part, so we take the real-part.
    # a real test function would suffice for a coeff <= 11.
    val = 11 + 72*np.sin(x)*np.cos(x)
    return np.sqrt(np.maximum(val,0))

def weight_fxn(i,x,t): #l_i(t)
    poly = 1
    for j in range(len(x)):
        if (i-1 == j):
            poly = poly*1
        else:
            poly = poly*((t - x[j])/(x[i-1]-x[j]))
    return poly

def get_weights_and_abscissas(n):
    roots = np.roots(spec.legendre(n))
    weights = []
    for i in range(1,n+1):
        wi = sp.quad(lambda t: weight_fxn(i,roots,t),-1,1)
        weights.append(wi[0])
    return weights, roots

# 5-point Gaussian Quadrature
weights5, roots5 = get_weights_and_abscissas(5)
#print(weights5, roots5) # optional

x = roots5
#print(sum(weights5*test_function(x))) # optional


def my_single_integral(f, A, B):
    weights, roots = get_weights_and_abscissas(5)
    roots = A + (B-A)/2 * (roots+1)
    return (B-A)/2 * sum(weights*test_function(roots))

# to test the 1-D Gaussian quadrature for [0,2].
#print(my_single_integral(test_function, 0, 2))


print("single integral test solution: ", my_single_integral(test_function, 0,np.pi))

def g(x): return 0
def h(x): return 1

def test_fxn2(x): return np.sin(5*x)

def my_double_integral(u, A, B, G, H):
    weights, roots = get_weights_and_abscissas(5)
    
    sum = 0
    
    x = A + (B-A)/2 * (roots + 1)
    
    for x_i,i in zip(x,range(5)):
        if not callable(H):
            d = H
        else:
            d = H(x_i)
        if not callable(G):
            c = G
        else:
            c = G(x_i)
        
        def F(y): return u(x_i, y)

        sum += weights[i]*my_single_integral(F, c, d)
        
    return (B-A)/2 * sum

print("double integral test solution: ", my_double_integral(test_fxn2, 0, 2, -1, 1))














