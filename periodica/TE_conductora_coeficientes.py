# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 17:56:24 2019

@author: Maria Sol Vidal
"""
# MODO TE perfectamente conductora 
import numpy as np
import cmath
#%% #defino parametros 
N=50
d=500 # en nanometros (0.5 micrometroa)
k=2*np.pi/d #k mayuscula 
a=645 #en nanometros
k_1=2*np.pi/a
tita=0 #incidencia normal 
alpha=k_1*np.sin(tita)
beta=k_1*np.cos(tita)
alpha_n = lambda x: alpha+x*k
def beta_n(x):
    return np.asarray(list(map(cmath.sqrt, k_1**2-alpha_n(x)**2)))
f = lambda x: 0 #perfil de la red 

#%%

aux = lambda x: k_1**2-alpha_n(x)**2
#%%
r=np.linspace(-N,N,2*N+1)
xp=np.linspace(0,d,2*N+1,endpoint=False) #preguntar esto 
#%%
def phi_p(x):
    return np.e**(1j*(alpha_n(r)*x+beta_n(r)*f(x)))
#%%
def ei(x):
    return -np.e**(1j*(alpha*x-beta*f(x)))
#%%
m=np.zeros((2*N+1,2*N+1),dtype=complex)
for i in range(len(xp)):
    m[i,:] = phi_p(xp[i])
#%%
b=ei(xp)
s=np.linalg.solve(m, b)
#%%
>>> lst = [-4, 3, -8, -9]
>>> map(cmath.sqrt, lst)
[2j, (1.7320508075688772+0j), 2.8284271247461903j, 3j]
>>> [cmath.sqrt(x) for x in lst]
[2j, (1.7320508075688772+0j), 2.8284271247461903j, 3j]