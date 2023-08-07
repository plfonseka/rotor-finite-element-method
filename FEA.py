# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 19:21:22 2020

@author: Lucas
"""

import numpy as np
import matplotlib.pyplot as plt
import vibration_toolbox as vbt
from plot_freq import plot_freq_resp
from mancal1 import mancal1
from mancal2 import mancal2

# Problem parameters
g = 9.8
E = 200E9
de = 10E-3
le = 800E-3
dd = 90E-3
ld = 15E-3
Iyy = (np.pi/64)*de**4
Izz = Iyy
a = 2/5*le
b = le-a
rho = 7850
vol_d = (np.pi*((dd/2)**2))*ld
m = rho*vol_d
P = m*g
I = np.pi*(de**4)/64
Id = (m/12)*((3/4)*dd**2 + ld**2)
#Id = m*((dd/2)**2)/4
Ip = (m*dd**2)/8
# Ip = m*((dd/2)**2)/2
A = np.pi*(de/2)**2

dm = 30E-3
rm = dm/2
lm = 20E-3
cr = 90E-6
neta = 0.07 
F1 = 0.6*P
F2 = 0.4*P

f4 = []
f5 = []
f6 = []
f7 = []
f8 = []
f9 = []
f10 = []
f11 = []
f12 = []
f13 = []
f14 = []
f15 = []
f16 = []
f17 = []
f18 = []
f19 = []
f20 = []
f21 = []
f22 = []
f23 = []
omega = []
omega2 = []

Nel = 5 # Number of mesh elements
NumNodes = Nel + 1 # Number of mesh nodes
ndof_node = 4
he = le/Nel # Element length
Ndofs = ndof_node * NumNodes # Total number of mesh dofs
X = np.linspace(0,le,NumNodes) # Nodal coordinates
x_aux = np.linspace(0,5,NumNodes)

# Matriz de Rigidez do elemento
Ke = (E*Iyy/(he**3)) * np.array([[12,0,0,6*he,-12,0,0,6*he],
                                 [0,12,-6*he,0,0,-12,-6*he,0],
                                 [0,-6*he,4*he**2,0,0,6*he,2*he**2,0],
                                 [6*he,0,0,4*he**2,-6*he,0,0,2*he**2],
                                 [-12,0,0,-6*he,12,0,0,-6*he],
                                 [0,-12,6*he,0,0,12,6*he,0],
                                 [0,-6*he,2*he**2,0,0,6*he,4*he**2,0],
                                 [6*he,0,0,2*he**2,-6*he,0,0,4*he**2]])

# Matriz de Massa do elemento
Me = ((rho*A*he)/420) * np.array([[156,0,0,22*he,54,0,0,-13*he],
                                  [0,156,-22*he,0,0,54,13*he,0],
                                  [0,-22*he,4*he**2,0,0,-13*he,-3*he**2,0],
                                  [22*he,0,0,4*he**2,13*he,0,0,-3*he**2],
                                  [54,0,0,13*he,156,0,0,-22*he],
                                  [0,54,-13*he**2,0,0,156,22*he,0],
                                  [0,13*he,-3*he**2,0,0,22*he,4*he**2,0],
                                  [-13*he,0,0,-3*he**2,-22*he,0,0,4*he**2]])

# Matriz Giroscópica do elemento
Ge = ((rho*A*de**2)/240*he) * np.array([[0,-36,3*he,0,0,36,3*he,0],
                                        [36,0,0,3*he,-36,0,0,3*he],
                                        [-3*he,0,0,-4*he**2,3*he,0,0,he**2],
                                        [0,-3*he,4*he**2,0,0,3*he,-he**2,0],
                                        [0,36,-3*he,0,0,-36,-3*he,0],
                                        [-36,0,0,-3*he,36,0,0,-3*he],
                                        [-3*he,0,0,he**2,3*he,0,0,-4*he**2],
                                        [0,-3*he,-he**2,0,0,3*he,4*he**2,0]])

Kg = np.zeros((Ndofs,Ndofs))
Mg = np.zeros((Ndofs,Ndofs))
Gg = np.zeros((Ndofs,Ndofs))
   
# Montagem das Matrizes Globais
for i in range(Nel):
    K_temp = np.zeros((Ndofs,Ndofs))
    M_temp = np.zeros((Ndofs,Ndofs))
    G_temp = np.zeros((Ndofs,Ndofs))
    K_temp[4*i:4*i+8, 4*i:4*i+8] = Ke
    M_temp[4*i:4*i+8, 4*i:4*i+8] = Me
    G_temp[4*i:4*i+8, 4*i:4*i+8] = Ge
    Kg += K_temp
    Mg += M_temp
    Gg += G_temp

alpha = 1E-4
Ce = alpha*Kg

'''
Adicionar efeitos do disco no sistema (nó 3)
'''

# massa do disco
Md = np.array([[m,0,0,0],
              [0,m,0,0],
              [0,0,Id,0],
              [0,0,0,Id]])

# efeito giroscópico
Gd = np.array([[0,0,0,0],
              [0,0,0,0],
              [0,0,0,-Ip],
              [0,0,Ip,0]])

Mg[8:12,8:12] += Md
Gg[8:12,8:12] += Gd

'''
Procedimento para plotagem do Diagrama de Campbell
'''
n = 4
i = 0
'''
  Adicionar efeitos dos Mancais no sistema (nós 1 e 6)
'''
mancal_1 = mancal1()

kyy1 = mancal_1[0]
kyz1 = mancal_1[1]
kzy1 = mancal_1[2]
kzz1 = mancal_1[3]

cyy1 = mancal_1[4]
cyz1 = mancal_1[5]
czy1 = mancal_1[6]
czz1 = mancal_1[7]

mancal_2 = mancal2()

kyy2 = mancal_2[0]
kyz2 = mancal_2[1]
kzy2 = mancal_2[2]
kzz2 = mancal_2[3]

cyy2 = mancal_2[4]
cyz2 = mancal_2[5]
czy2 = mancal_2[6]
czz2 = mancal_2[7]    

for W in np.arange(1,2000,5):

    # Matriz de rigidez p/ o mancal 1 (nó 1)
    Km1 = np.array([[kyy1[i],kyz1[i],0,0],                  
                    [kzy1[i],kzz1[i],0,0],
                    [0,0,0,0],
                    [0,0,0,0]])
    # Matriz de rigidez p/ o mancal 2 (nó 6)
    Km2 = np.array([[kyy2[i],kyz2[i],0,0],
                    [kzy2[i],kzz2[i],0,0],
                    [0,0,0,0],
                    [0,0,0,0]])
    
    # Matriz de amortecimento p/ o mancal 1 (nó 1)
    Cm1 = np.array([[cyy1[i],cyz1[i],0,0],
                    [czy1[i],czz1[i],0,0],
                    [0,0,0,0],
                    [0,0,0,0]])
    # Matriz de amortecimento p/ o mancal 2 (nó 6)
    Cm2 = np.array([[cyy2[i],cyz2[i],0,0],
                    [czy2[i],czz2[i],0,0],
                    [0,0,0,0],
                    [0,0,0,0]])

    # atualização dos mancais na Matriz de Rigidez Global
    Kg[0:n,0:n] += Km1
    Kg[20:20+n,20:20+n] += Km2
    
    # atualização dos mancais na Matriz de Amortecimento Global
    Ce[0:n,0:n] += Cm1
    Ce[20:20+n,20:20+n] += Cm2
    
    Cg = Ce + W*Gg
    
    sys = vbt.VibeSystem(Mg,Cg,Kg)  
    freq = sys.wd*0.16
    
    f4.append(freq[4])
    f5.append(freq[5])
    # f6.append(freq[6])
    # f7.append(freq[7])
    f8.append(freq[8])
    f9.append(freq[9])
    f10.append(freq[10])
    f11.append(freq[11])
    f12.append(freq[12])
    f13.append(freq[13])
    f14.append(freq[14])
    f15.append(freq[15])
    # # f16.append(freq[16])
    # # f17.append(freq[17])
    # # f18.append(freq[18])
    # # f19.append(freq[19])
    # # f20.append(freq[20])
    # # f21.append(freq[21])
    # # f22.append(freq[22])
    # # f23.append(freq[23])
    
    omega.append(W*0.16)
    omega2.append((2*W)*0.16)
    i = i+1

plot_freq_resp(sys)

fig1, ax1 = plt.subplots()
ax1.set_title('Diagrama de Campbell - Sistema Amortecido')
ax1.set_xlabel('Rotação [Hz]')
ax1.set_ylabel('Frequência Natural Amortecida ($\omega_d$)')

ax1.plot(omega,f4,'k')
ax1.plot(omega,f5,'grey')
# ax1.plot(omega,f6)
# ax1.plot(omega,f7)
ax1.plot(omega,f8,'k')
ax1.plot(omega,f9,'grey')
ax1.plot(omega,f10,'k')
ax1.plot(omega,f11,'grey')
ax1.plot(omega,f12,'k')
ax1.plot(omega,f13,'grey')
ax1.plot(omega,f14,'k')
ax1.plot(omega,f15,'grey')
# ax1.plot(omega,f16)
# ax1.plot(omega,f17)
# ax1.plot(omega,f18)
# ax1.plot(omega,f19)
# ax1.plot(omega,f20)
# ax1.plot(omega,f21)
# ax1.plot(omega,f22)
# ax1.plot(omega,f23)

ax1.plot(omega,omega,'b--')
ax1.plot(omega,omega2,'r--')
ax1.margins(x=0.0)
# plt.axis('image')
# plt.savefig('frf_disco.png',dpi=600,bbox_inches = 'tight')
