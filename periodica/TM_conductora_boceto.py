# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 17:56:24 2019

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
alpha_n = lambda x: alpha+x*k
betasca_n = lambda x: cmath.sqrt(k_1**2-alpha_n(x)**2)
def beta_n(x):
    return np.array(list(map(cmath.sqrt, k_1**2-alpha_n(x)**2)))
f = lambda x: 0 #perfil de la red 

#%%

def ab(x,t):
    return np.array([[alpha_n(i),betasca_n(i)]*(np.e**(1j*(alpha_n(i)*t+betasca_n(i)*f(t))))  for i in x])
#%%
def ab2(r,t):
    h=np.zeros((2*N+1,2),dtype=complex)
    for i in range(2*N+1):
        h[i,:]=np.array([alpha_n(r[i]),betasca_n(r[i])])*np.e**(1j*(alpha_n(r[i])*t+betasca_n(r[i])*f(t)))
    return h 
#%%
    
m=ab2(r,np.pi)
b=normal(np.cos,np.pi)
j=m.dot(b)
    #%%
j=np.array([[alpha_n(i),betasca_n(i)]*np.e**(1j*(alpha_n(i)*t+betasca_n(i)*f(t)))  for i in x])
#%%
def normal(dy,t):
    return ((1+dy(t)**2)**(-1/2))*np.array([-dy(t),1])
#%%
vn=np.zeros((2*N+1,1),dtype=complex)
for i in range(len(xp)):
    m[i,:] = phi_p(xp[i])




#%%
def pi(dy,t,n):
    return np.array([np.dot(normal(dy,t),i)for i in ab(n,t)])
#%%
cuadratica = lambda x: x**2  
#%%
r=np.linspace(-N,N,2*N+1)
xp=np.linspace(0,d,2*N+1,endpoint=False) #preguntar esto 
#%%
np.dot(normal(cuadratica,0),normal(cuadratica,1))
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
#%%
>>> lst = [-4, 3, -8, -9]
>>> map(cmath.sqrt, lst)
[2j, (1.7320508075688772+0j), 2.8284271247461903j, 3j]
>>> [cmath.sqrt(x) for x in lst]
[2j, (1.7320508075688772+0j), 2.8284271247461903j, 3j]