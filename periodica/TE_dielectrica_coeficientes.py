# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 19:09:28 2019

@author: Maria Sol Vidal
"""
import numpy as np
import cmath
#%% #defino parametros 
N=50
nu=1.45 # indice del medio de abajo 
d=500 # en nanometros (0.5 micrometroa)
k=2*np.pi/d #k mayuscula 
a=645 #en nanometros #longitud de onda incidente 
k_1=2*np.pi/a #ene la region incidente hay aire
tita=np.pi/4 #incidencia normal  #angulo de incidencia
alpha=k_1*np.sin(tita)
beta=k_1*np.cos(tita)
#%% #Funciones
r=np.linspace(-N,N,2*N+1)
alpha_n = lambda n: alpha+n*k
betasca_n = lambda n: cmath.sqrt(k_1**2-alpha_n(n)**2)
betasca_n2 = lambda n: cmath.sqrt((k_1*nu)**2-alpha_n(n)**2)
def beta_n(n):
    return np.array(list(map(cmath.sqrt, k_1**2-alpha_n(n)**2)))
def beta_n2(n):
    return np.array(list(map(cmath.sqrt, (k_1*nu)**2-alpha_n(n)**2)))
f = lambda x: 0 #perfil de la red 
def phi_p(x):
    return np.e**(1j*(alpha_n(r)*x+beta_n(r)*f(x)))
def phi_p2(x):
    return np.e**(1j*(alpha_n(r)*x-beta_n2(r)*f(x)))
def ei1(x):
    return -np.e**(1j*(alpha*x-beta*f(x)))
def ab2(r,x):
    h=np.zeros((len(r),2),dtype=complex)
    for i in range(len(r)):
        h[i,:]=np.array([alpha_n(r[i]),betasca_n(r[i])])*np.e**(1j*(alpha_n(r[i])*x+betasca_n(r[i])*f(x)))
    return h 
def ab3(r,x):
    h=np.zeros((len(r),2),dtype=complex)
    for i in range(len(r)):
        h[i,:]=np.array([alpha_n(r[i]),-betasca_n2(r[i])])*np.e**(1j*(alpha_n(r[i])*x-betasca_n2(r[i])*f(x)))
    return h 
def normal(df,x):
    return ((1+df(x)**2)**(-1/2))*np.array([-df(x),1])
def ei2(x):
    return -np.dot(normal(f,x),np.array([alpha,-beta]))*np.e**(1j*(alpha*x-beta*f(x)))
#%%
def coef_TE(N,d):
    r=np.linspace(-N,N,2*N+1) #vector con Ns 
    xp=np.linspace(0+10**(-6),d-10**(-6),2*N+1) #preguntar esto
    b=np.zeros((4*N+2),dtype=complex)
    b[0:2*N+1]=ei1(xp)
    for i in range(len(xp)):
        b[i+2*N+1] = ei2(xp[i])
    m=np.zeros((4*N+2,4*N+2),dtype=complex)
#lleno filas 0 a 2N y columnas 0 a 2N
    for i in range(len(xp)):
        m[i,0:2*N+1] = phi_p(xp[i])
#lleno filas 0 a 2N y columnas 2N+1, 4N+1
    for i in range(len(xp)):
        m[i,2*N+1:4*N+2] = -phi_p2(xp[i])
    for i in range(len(r)):
        m[i+2*N+1,0:2*N+1] = ab2(r,xp[i]).dot(normal(f,xp[i]))
    for i in range(len(r)):
        m[i+2*N+1,2*N+1:4*N+2] = -ab3(r,xp[i]).dot(normal(f,xp[i]))
    return np.linalg.solve(m,b)
#%%
np.conj(1j)
#%%
coef_TE(N,d)




