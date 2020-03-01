

import numpy as np
import cmath
#%%
import matplotlib.pyplot as plt
#%% #defino parametros 
N=100
nu=2.25**(1/2) # indice del medio de abajo 
d=2*np.pi # en nanometros (0.5 micrometroa)
a=0.5*np.pi#en nanometros #longitud de onda incidente 
tita=0 #incidencia normal  #angulo de incidencia
#tita=np.arange(0,np.pi/2,0.01)
f = lambda x: np.cos(2*np.pi*x/a)*0.005*a/2 #perfil de la red 
df = lambda x: -np.sin(2*np.pi*x/a)*np.pi*0.005
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
    return h 
def normal(df,x):
    return ((1+df(x)**2)**(-1/2))*np.array([-df(x),1])
def ei2(f,df,x,alpha,beta):
    return -np.dot(normal(df,x),np.array([alpha,-beta]))*np.e**(1j*(alpha*x-beta*f(x)))
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
    efici2=(h2[wh22]*np.absolute(x[wh2])**2)/(k_1**2-alpha**2)**(1/2)
    return [efici,efici2,orden_refle(wh,N),orden_trans(wh2,N)]

#%%
def orden_efi(efi,orden):
    dict= {efi[i] : orden[i] for i in range(len(efi))}
    return dict

#%%
    #nu*sen(tita1n)=sen(tita)+nlambda/d
def ecuacion(nu,tita,n,lambdaa,d):
    return (np.sin(tita)+n*lambdaa/d)/nu
#%%
lista=[]
chequeo=[]
for ang in tita:
    efi_refle, efi_trans, orden_refl, orden_tran=efi_TE(N,nu,d,a,ang,f,df)
    dic=orden_efi(efi_refle,orden_refl)
    chequeo.append(sum(efi_refle)+sum(efi_trans))
    for i in efi_refle:
        if dic[i]==-1:
            lista.append([i,ang])
            #%%
angulo=[]
efi_menos1=[]
for i in lista:
    angulo.append(i[1])
    efi_menos1.append(i[0])
    #%%
#%%
plt.plot(angulo,efi_menos1,'go',markersize=3,label='Método de Rayleigh')
plt.title("Eficiencia para el orden m=-1 para una interfaz sinusoidal.nu=1.45, lambda=645nm ")
plt.ylabel("Eficiencia")
plt.xlabel("Ángulo de incidencia (radianes)")
plt.legend()
#%%
che=np.array(chequeo).real
plt.plot(tita,che,'go',markersize=3,label="Balance energético")
plt.ticklabel_format(useOffset=False,style='sci')
plt.ylim((0.9999999999998, 1.0000000000002))
plt.title("Suma de eficiencia para una interfaz sinusoidal.nu=1.45, lambda=645nm ")
#plt.ylabel("Eficiencia")
plt.xlabel("Ángulo de incidencia (radianes)")
plt.legend()
#%%
N=50
nu=1.45 # indice del medio de abajo 
d=200*np.pi # en nanometros (0.5 micrometroa)
a=645 #en nanometros #longitud de onda incidente 
tita=np.pi/4 #incidencia normal  #angulo de incidencia
f = lambda x: np.sin(x/100) #perfil de la red 
df = lambda x: (1/100)*np.cos(x/100)
long=np.linspace(400,700,50)
lista2=[]
chequeo2=[]
for lo in long:
    efi_refle, efi_trans, orden_refl, orden_tran=efi_TE(N,nu,d,lo,tita,f,df)
    dic=orden_efi(efi_refle,orden_refl)
    chequeo2.append(sum(efi_refle)+sum(efi_trans))
    for i in efi_refle:
        if dic[i]==-1:
            lista2.append([i,lo])
#%%
lam=[]
efi_menos1=[]
for i in lista2:
    lam.append(i[1])
    efi_menos1.append(i[0])          
            
#%%
plt.plot(lam,efi_menos1,'go',markersize=3,label='Método de Rayleigh')
plt.title("Eficiencia para el orden m=-1 para una interfaz sinusoidal.nu=1.45, tita=pi/4 ")
plt.ylabel("Eficiencia")
plt.xlabel("Longitud de onda (nm)")
plt.legend()       
            
            
            
            
#%%

efi_refle, efi_trans, orden_refl, orden_tran=efi_TE(N,nu,d,a,tita[0],f,df)
dic=orden_efi(efi_refle,orden_refl)
   # chequeo.append(sum(efi_refle)+sum(efi_trans))
for i in efi_refle:
    if dic[i]==-1:
        lista.append([i,ang])      