# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:47:06 2019

@author: Maria Sol Vidal
"""

import numpy as np
import cmath
#%% #defino parametros 
N=100
nu=2.25**(1/2) # indice del medio de abajo 
d=2*np.pi # en nanometros (0.5 micrometroa)
a=d/1.7 #en nanometros #longitud de onda incidente 
tita=30*np.pi/180 #incidencia normal  #angulo de incidencia
#tita=np.arange(0,np.pi/2,0.01)
f = lambda x: 0.1*d*np.cos(x)/2 #perfil de la red 
df = lambda x: -0.1*d*np.sin(x)/2
#f = lambda x: 0 #perfil de la red
#long=np.linspace(400,700,60)


#%%
alpha_n = lambda alpha,n,k: alpha+n*k
betasca_n = lambda alpha,n,k,k_1: cmath.sqrt(k_1**2-alpha_n(alpha,n,k)**2)
betasca_n2 = lambda alpha,n,k,k_1,nu: cmath.sqrt((k_1*nu)**2-alpha_n(alpha,n,k)**2)
def beta_n(alpha,n,k,k_1):
    return np.array(list(map(cmath.sqrt, k_1**2-alpha_n(alpha,n,k)**2)))
def beta_n2(alpha,n,k,k_1,nu):
    return np.array(list(map(cmath.sqrt, (k_1*nu)**2-alpha_n(alpha,n,k)**2)))
def beta_np(alpha,n,k,k_1):
    return k_1**2-alpha_n(alpha,n,k)**2
def beta_n2p(alpha,n,k,k_1,nu):
    return (k_1*nu)**2-alpha_n(alpha,n,k)**2
def phi_p(alpha,r,x,k,f,k_1):
    return np.e**(1j*(alpha_n(alpha,r,k)*x+beta_n(alpha,r,k,k_1)*f(x)))
def phi_p2(alpha,r,x,k,f,k_1,nu):
    return np.e**(1j*(alpha_n(alpha,r,k)*x-beta_n2(alpha,r,k,k_1,nu)*f(x)))
def ei1(x,alpha,beta,f):
    return -np.e**(1j*(alpha*x-beta*f(x)))
def ab2(alpha,r,x,k,f,k_1):
    h=np.zeros((len(r),2),dtype=complex)
    for i in range(len(r)):
        h[i,:]=np.array([alpha_n(alpha,r[i],k),betasca_n(alpha,r[i],k,k_1)])*np.e**(1j*(alpha_n(alpha,r[i],k)*x+betasca_n(alpha,r[i],k,k_1)*f(x)))
    return h 
def ab3(alpha,r,x,k,f,k_1,nu):
    h=np.zeros((len(r),2),dtype=complex)
    for i in range(len(r)):
        h[i,:]=np.array([alpha_n(alpha,r[i],k),-betasca_n2(alpha,r[i],k,k_1,nu)])*np.e**(1j*(alpha_n(alpha,r[i],k)*x-betasca_n2(alpha,r[i],k,k_1,nu)*f(x)))
    return h/nu**2 
def normal(df,x):
    return ((1+df(x)**2)**(-1/2))*np.array([-df(x),1])
def ei2(f,df,x,alpha,beta):
    return -np.dot(normal(df,x),np.array([alpha,-beta]))*np.e**(1j*(alpha*x-beta*f(x)))
#tengo que cambiar ,ab3 porq e=1 para arriba 
#%%
orden_refle = lambda indice,N: indice-N
orden_trans = lambda indice,N: indice-(3*N+1)
#%%
def efi_TE(N,nu,d,a,tita,f,df):
    k=2*np.pi/d #k mayuscula 
    k_1=2*np.pi/a #ene la region incidente hay aire
    alpha=k_1*np.sin(tita)
    beta=k_1*np.cos(tita)
    r=np.linspace(-N,N,2*N+1) #vector con Ns 
   # xp=np.linspace(0+10**(-6),d-10**(-6),2*N+1) #preguntar esto
    xp=np.linspace(-d/2+10**(-6),d/2-10**(-6),2*N+1) #preguntar esto
    b=np.zeros((4*N+2),dtype=complex)
    b[0:2*N+1]=ei1(xp,alpha,beta,f)
    for i in range(len(xp)):
        b[i+2*N+1] = ei2(f,df,xp[i],alpha,beta)
    m=np.zeros((4*N+2,4*N+2),dtype=complex)
#lleno filas 0 a 2N y columnas 0 a 2N
    for i in range(len(xp)):
        m[i,0:2*N+1] = phi_p(alpha,r,xp[i],k,f,k_1)
#lleno filas 0 a 2N y columnas 2N+1, 4N+1
    for i in range(len(xp)):
        m[i,2*N+1:4*N+2] = -phi_p2(alpha,r,xp[i],k,f,k_1,nu)
    for i in range(len(r)):
        m[i+2*N+1,0:2*N+1] = ab2(alpha,r,xp[i],k,f,k_1).dot(normal(df,xp[i]))
    for i in range(len(r)):
        m[i+2*N+1,2*N+1:4*N+2] = -ab3(alpha,r,xp[i],k,f,k_1,nu).dot(normal(df,xp[i]))
    x=np.linalg.solve(m,b)
    p1=beta_np(alpha,r,k,k_1)
    p2=beta_n2p(alpha,r,k,k_1,nu)
    h=beta_n(alpha,r,k,k_1)
    h2=beta_n2(alpha,r,k,k_1,nu)
    wh=np.where(p1>0)[0]
    wh2=np.where(p2>0)[0]+2*N+1
    wh22=np.where(p2>0)[0]
    efici=(h[wh]*np.absolute(x[wh])**2)/(k_1**2-alpha**2)**(1/2)
    efici2=(h2[wh22]*np.absolute(x[wh2])**2)/(((k_1**2-alpha**2)**(1/2))*(nu**2))
    return [efici,efici2,orden_refle(wh,N),orden_trans(wh2,N)]

#%%
