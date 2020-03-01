# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 12:16:24 2019

@author: Maria Sol Vidal
"""
import numpy as np
import cmath

#%% #defino parametros 
N=50
d=2000000*np.pi # en nanometros (0.5 micrometroa)
a=0.8*d #en nanometros #longitud de onda incidente 
tita=0 #incidencia normal  #angulo de incidencia
#tita=np.arange(0,np.pi/2,0.01)
f = lambda x: 0.2*d*np.sin(x/1000000) #perfil de la red 
#df = lambda x: 0.4*d*np.cos(x/100)/100
#f = lambda x: 0 #perfil de la red
#%%
alpha_n = lambda alpha,n,k: alpha+n*k
def beta_n(alpha,n,k,k_1):
    return np.array(list(map(cmath.sqrt, k_1**2-alpha_n(alpha,n,k)**2)))
def phi_p(alpha,r,x,k,f,k_1):
    return np.e**(1j*(alpha_n(alpha,r,k)*x+beta_n(alpha,r,k,k_1)*f(x)))
def ei(x,alpha,beta,f):
    return -np.e**(1j*(alpha*x-beta*f(x)))
def beta_np(alpha,n,k,k_1):
    return k_1**2-alpha_n(alpha,n,k)**2
#%%
orden_refle = lambda indice,N: indice-N


#%%
def efi_TE(N,d,a,tita,f):
    k=2*np.pi/d #k mayuscula 
    k_1=2*np.pi/a #ene la region incidente hay aire
    alpha=k_1*np.sin(tita)
    beta=k_1*np.cos(tita)
    r=np.linspace(-N,N,2*N+1) #vector con Ns 
    xp=np.linspace(0+10**(-6),d-10**(-6),2*N+1) #preguntar esto
    #xp=np.linspace(-d/2+10**(-6),d/2-10**(-6),2*N+1) #preguntar esto
    b=ei(xp,alpha,beta,f)
    m=np.zeros((2*N+1,2*N+1),dtype=complex)
    for i in range(len(xp)):
        m[i,:] = phi_p(alpha,r,xp[i],k,f,k_1)
    x=np.linalg.solve(m,b)
    p1=beta_np(alpha,r,k,k_1)
    h=beta_n(alpha,r,k,k_1)
    wh=np.where(p1>0)[0]
    efici=(h[wh]*np.absolute(x[wh])**2)/(k_1**2-alpha**2)**(1/2)
    return [efici,orden_refle(wh,N)]

#%%
resul=efi_TE(N,d,a,tita,f)
#%%
k=2*np.pi/d #k mayuscula 
k_1=2*np.pi/a #ene la region incidente hay aire
alpha=k_1*np.sin(tita)
beta=k_1*np.cos(tita)
r=np.linspace(-N,N,2*N+1) #vector con Ns 
xp=np.linspace(0+10**(-6),d-10**(-6),2*N+1) #preguntar esto
#%%
b=ei(xp,alpha,beta,f)
m=np.zeros((2*N+1,2*N+1),dtype=complex)
for i in range(len(xp)):
    m[i,:] = phi_p(alpha,r,xp[i],k,f,k_1)
x=np.linalg.solve(m,b)
#%%
p1=beta_np(alpha,r,k,k_1)
#%%
h=beta_n(alpha,r,k,k_1)
wh=np.where(p1>0)[0]
#%%
efici=(h[wh]*np.absolute(x[wh])**2)/(k_1**2-alpha**2)**(1/2)

