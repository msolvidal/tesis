# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 13:40:31 2019

@author: Maria Sol Vidal
"""

# En mis programas ni=1 siempre  y mu=1
#np.cos() esta en radianes
#%%
import numpy as np
#%%
def costitat(ni,nt,titai):
    return (1-(ni/nt*np.sin(titai))**2)**(1/2)
#%% #TE
def rper(ni,nt,titai):
    return (ni*np.cos(titai)-nt*costitat(ni,nt,titai))/(ni*np.cos(titai)+nt*costitat(ni,nt,titai))
#%% #TM
def r_paralel(ni,nt,titai):
    return (nt*np.cos(titai)-ni*costitat(ni,nt,titai))/(nt*np.cos(titai)+ni*costitat(ni,nt,titai))
#%% #TM
def t_paralel(ni,nt,titai):
    return (2*ni*np.cos(titai))/(nt*np.cos(titai)+ni*costitat(ni,nt,titai))
#%% #TE
def tper(ni,nt,titai):
    return (2*ni*np.cos(titai))/(ni*np.cos(titai)+nt*costitat(ni,nt,titai)) 
#%%   
lista2=[]
for ang in tita:
    lista2.append([rper(1,nu,ang),tper(1,nu,ang)])
#%%
tita=np.pi/4 
val_fres=[rper(1,nu,tita),tper(1,nu,tita)]