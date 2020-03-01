# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:44:22 2019

@author: Maria Sol Vidal
"""

# MODO TM perfectamente conductora 
import numpy as np
import cmath
#%% #defino parametros 
N=50
d=500 # en nanometros (0.5 micrometroa)
k=2*np.pi/d #k mayuscula 
a=645 #en nanometros #longitud de onda incidente 
k_1=2*np.pi/a #ene la region incidente hay aire
tita=0 #incidencia normal  #angulo de incidencia
alpha=k_1*np.sin(tita)
beta=k_1*np.cos(tita)
alpha_n = lambda n: alpha+n*k
betasca_n = lambda n: cmath.sqrt(k_1**2-alpha_n(n)**2)
def beta_n(n):
    return np.array(list(map(cmath.sqrt, k_1**2-alpha_n(n)**2)))
f = lambda x: 0 #perfil de la red 
#%%
def ab2(r,x):
    h=np.zeros((len(r),2),dtype=complex)
    for i in range(len(r)):
        h[i,:]=np.array([alpha_n(r[i]),betasca_n(r[i])])*np.e**(1j*(alpha_n(r[i])*x+betasca_n(r[i])*f(x)))
    return h 
#%%
def normal(df,x):
    return ((1+df(x)**2)**(-1/2))*np.array([-df(x),1])
#%%
r=np.linspace(-N,N,2*N+1)
xp=np.linspace(0+10**(-6),d-10**(-6),2*N+1) #preguntar esto 
#%%
m=np.zeros((2*N+1,2*N+1),dtype=complex)
for i in range(len(r)):
    m[i,:] = ab2(r,xp[i]).dot(normal(f,xp[i]))
    
#%%
def ei(x):
    return -np.dot(normal(f,x),np.array([alpha,-beta]))*np.e**(1j*(alpha*x-beta*f(x)))
#%%
b=np.zeros((2*N+1,1),dtype=complex)
for i in range(len(xp)):
    b[i,:] = ei(xp[i])
#%%
s=np.linalg.solve(m,b)
    
    
    
    