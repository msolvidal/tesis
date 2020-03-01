# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 20:37:09 2019

@author: Maria Sol Vidal
"""

import numpy as np
import cmath
#%% #defino parametros 
N=50
nu=1.45 # indice del medio de abajo 
d=200*np.pi # en nanometros (0.5 micrometroa)
a=645 #en nanometros #longitud de onda incidente 
#tita=np.pi/4 #incidencia normal  #angulo de incidencia
tita=np.arange(0,np.pi/2,0.01)
f = lambda x: np.sin(x/100) #perfil de la red 
df = lambda x: (1/100)*np.cos(x/100)
#f = lambda x: 0 #perfil de la red
#%%

alpha_n = lambda alpha,n,k: alpha+n*k
betasca_n = lambda alpha,n,k,k_1: cmath.sqrt(k_1**2-alpha_n(alpha,n,k)**2)
betasca_n2 = lambda alpha,n,k,k_1,nu: cmath.sqrt((k_1*nu)**2-alpha_n(alpha,n,k)**2)
def beta_n(alpha,n,k,k_1):
    return np.array(list(map(cmath.sqrt, k_1**2-alpha_n(alpha,n,k)**2)))
def beta_n2(alpha,n,k,k_1,nu):
    return np.array(list(map(cmath.sqrt, (k_1*nu)**2-alpha_n(alpha,n,k)**2)))
f = lambda x: 0 #perfil de la red 
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
    return h 
def normal(df,x):
    return ((1+df(x)**2)**(-1/2))*np.array([-df(x),1])
def ei2(f,x,alpha,beta):
    return -np.dot(normal(f,x),np.array([alpha,-beta]))*np.e**(1j*(alpha*x-beta*f(x)))
#%%
def coef_TE(N,nu,d,a,tita,f,df):
    k=2*np.pi/d #k mayuscula 
    k_1=2*np.pi/a #ene la region incidente hay aire
    alpha=k_1*np.sin(tita)
    beta=k_1*np.cos(tita)
    r=np.linspace(-N,N,2*N+1) #vector con Ns 
    xp=np.linspace(0+10**(-6),d-10**(-6),2*N+1) #preguntar esto
    b=np.zeros((4*N+2),dtype=complex)
    b[0:2*N+1]=ei1(xp,alpha,beta,f)
    for i in range(len(xp)):
        b[i+2*N+1] = ei2(f,xp[i],alpha,beta)
    m=np.zeros((4*N+2,4*N+2),dtype=complex)
#lleno filas 0 a 2N y columnas 0 a 2N
    for i in range(len(xp)):
        m[i,0:2*N+1] = phi_p(alpha,r,xp[i],k,f,k_1)
#lleno filas 0 a 2N y columnas 2N+1, 4N+1
    for i in range(len(xp)):
        m[i,2*N+1:4*N+2] = -phi_p2(alpha,r,xp[i],k,f,k_1,nu)
    for i in range(len(r)):
        m[i+2*N+1,0:2*N+1] = ab2(alpha,r,xp[i],k,f,k_1).dot(normal(f,xp[i]))
    for i in range(len(r)):
        m[i+2*N+1,2*N+1:4*N+2] = -ab3(alpha,r,xp[i],k,f,k_1,nu).dot(normal(f,xp[i]))
    return np.linalg.solve(m,b) 

#%%
lista=[]
for ang in tita:
    g=coef_TE(N,nu,d,a,ang,f)
    lista.append([g[50],g[151]])
    #%%
a=np.asarray(lista)
a=a.real
#%%
b=np.asarray(lista2)
b=b.real
#%%
refl_ray=[]
refl_fres=[]
for i in range(len(tita)):
    refl_ray.append(a[i][0])
    refl_fres.append(b[i][0])
    #%%
import matplotlib.pyplot as plt
#%%
plt.plot(tita,refl_ray,'g',label='Método de Rayleigh')
plt.plot(tita,refl_fres,"r--",label='Coeficientes de Fresnel')
plt.title("Coeficiente R para una interfaz plana.nu=1.45, lambda=645nm ")
plt.ylabel("R")
plt.xlabel("'Ángulo de incidencia (radianes)")
plt.legend()
#%%
resta=np.asarray(refl_ray)-np.asarray(refl_fres)
#%%
max(resta)/0.6489857736444861
#%%
np.where(resta == np.amax(resta))
#%%
plt.plot(tita,resta,'go')
#plt.title("Resta de R por Fresnel y Rayleigh.Interfaz plana.nu=1.45, lambda=645nm ",loc='center')
plt.ylabel("Resta R por Rayleigh y Fresnel")
plt.xlabel("'Ángulo de incidencia (radianes)")
#%%
long=np.linspace(400,700,60)
lista_lon=[]
for lon in long:
    g=coef_TE(N,nu,d,lon,tita,f)
    lista_lon.append([g[50],g[151]])
    #%%
a_l=np.asarray(lista_lon)
a_l=a_l.real
#%%
refl_ray_l=[]
for i in range(len(long)):
    refl_ray_l.append(a_l[i][0])
    #%%

plt.plot(long,refl_ray_l,"bo",markersize=3,label='Método de Rayleigh')
plt.ticklabel_format(useOffset=False,style='sci')
plt.ylim((-0.28321870385778, -0.28321870385774))
#ax.ticklabel_format(useOffset=False, style='sci')
#plt.ticklabel_format(style='sci')
#plt.ticklabel_format(useOffset=False)
plt.axhline(y=val_fres[0], color='r',label='Coeficiente de Fresnel')
plt.title("Coeficiente R para una interfaz plana.nu=1.45, tita=pi/4 ")
plt.ylabel("R")
plt.xlabel("Longitud de onda (nm)")
plt.legend()
#%%
trans_ray_l=[]
for i in range(len(long)):
    trans_ray_l.append(a_l[i][1])
    #%%
plt.plot(long,trans_ray_l,"bo",markersize=3,label='Método de Rayleigh')
plt.ticklabel_format(useOffset=False,style='sci')
plt.ylim((0.7167812961422301, 0.7167812961422501))
#ax.ticklabel_format(useOffset=False, style='sci')
#plt.ticklabel_format(style='sci')
#plt.ticklabel_format(useOffset=False)
plt.axhline(y=val_fres[1], color='r',label='Coeficiente de Fresnel')
plt.title("Coeficiente T para una interfaz plana.nu=1.45, tita=pi/4 ")
plt.ylabel("T")
plt.xlabel("Longitud de onda (nm)")
plt.legend()
    

#%%
np.std(refl_ray_l)
#%%
np.mean(trans_ray_l)
#%%
max(np.asarray(refl_ray_l)-val_fres[0])/val_fres[0]
#%%
max(np.asarray(trans_ray_l)-val_fres[1])/val_fres[1]
#%%
max(np.asarray(trans_ray_l)-val_fres[1])

#%%
trans_ray=[]
trans_fres=[]
for i in range(len(tita)):
    trans_ray.append(a[i][1])
    trans_fres.append(b[i][1])
#%%
plt.plot(tita,trans_ray,'g',label='Método de Rayleigh')
plt.plot(tita,trans_fres,"r--",label='Coeficientes de Fresnel')
plt.title("Coeficiente T para una interfaz plana.nu=1.45, lambda=645nm ")
plt.ylabel("T")
plt.xlabel("Ángulo de incidencia (radianes)")
plt.legend()
#%%
resta=np.asarray(trans_ray)-np.asarray(trans_fres)
#%%
max(resta)/trans_fres[134]
#%%
np.where(resta == np.amax(resta))
#%%
r=abs(np.asarray(refl_ray))
t=np.asarray(trans_ray)
tot=r+t
#%%
plt.plot(tita,tot,'g')
plt.ticklabel_format(useOffset=False,style='sci')
plt.ylim((0.99999999999998, 1.000000000000025))
plt.title("Suma de R y T para una interfaz plana.nu=1.45, lambda=645nm ")
#plt.ylabel("T")
plt.xlabel("Ángulo de incidencia (radianes)")
plt.legend()

