# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 23:38:44 2020

@author: Lucas
"""

import numpy as np
import matplotlib.pyplot as plt
import vibration_toolbox as vbt
from mpl_toolkits import mplot3d
# from plot_freq import plot_freq_resp
from mancal1_temp import mancal1
from mancal2_temp import mancal2


# E = 200E9
# de = 10E-3
le = 800E-3
# dd = 90E-3
# ld = 15E-3
# Iyy = (np.pi/64)*de**4

# rho = 7850
# vol_d = (np.pi*((dd/2)**2))*ld
# m = rho*vol_d

# Id = (m/12)*((3/4)*dd**2 + ld**2)
# Ip = (m*dd**2)/8
# A = np.pi*(de/2)**2

# W = 100

Nel = 5 # Number of mesh elements
NumNodes = Nel + 1 # Number of mesh nodes
ndof_node = 4
he = le/Nel # Element length
Ndofs = ndof_node * NumNodes # Total number of mesh dofs

# # Matriz de Rigidez do elemento
# Ke = (E*Iyy/(he**3)) * np.array([[12,0,0,6*he,-12,0,0,6*he],
#                                  [0,12,-6*he,0,0,-12,-6*he,0],
#                                  [0,-6*he,4*he**2,0,0,6*he,2*he**2,0],
#                                  [6*he,0,0,4*he**2,-6*he,0,0,2*he**2],
#                                  [-12,0,0,-6*he,12,0,0,-6*he],
#                                  [0,-12,6*he,0,0,12,6*he,0],
#                                  [0,-6*he,2*he**2,0,0,6*he,4*he**2,0],
#                                  [6*he,0,0,2*he**2,-6*he,0,0,4*he**2]])

# # Matriz de Massa do elemento
# Me = ((rho*A*he)/420) * np.array([[156,0,0,22*he,54,0,0,-13*he],
#                                   [0,156,-22*he,0,0,54,13*he,0],
#                                   [0,-22*he,4*he**2,0,0,-13*he,-3*he**2,0],
#                                   [22*he,0,0,4*he**2,13*he,0,0,-3*he**2],
#                                   [54,0,0,13*he,156,0,0,-22*he],
#                                   [0,54,-13*he**2,0,0,156,22*he,0],
#                                   [0,13*he,-3*he**2,0,0,22*he,4*he**2,0],
#                                   [-13*he,0,0,-3*he**2,-22*he,0,0,4*he**2]])

# # Matriz Giroscópica do elemento
# Ge = ((rho*A*de**2)/240*he) * np.array([[0,-36,3*he,0,0,36,3*he,0],
#                                         [36,0,0,3*he,-36,0,0,3*he],
#                                         [-3*he,0,0,-4*he**2,3*he,0,0,he**2],
#                                         [0,-3*he,4*he**2,0,0,3*he,-he**2,0],
#                                         [0,36,-3*he,0,0,-36,-3*he,0],
#                                         [-36,0,0,-3*he,36,0,0,-3*he],
#                                         [-3*he,0,0,he**2,3*he,0,0,-4*he**2],
#                                         [0,-3*he,-he**2,0,0,3*he,4*he**2,0]])

# Kg = np.zeros((Ndofs,Ndofs))
# Mg = np.zeros((Ndofs,Ndofs))
# Gg = np.zeros((Ndofs,Ndofs))
   
# K_temp = np.zeros((Ndofs,Ndofs))
# M_temp = np.zeros((Ndofs,Ndofs))
# G_temp = np.zeros((Ndofs,Ndofs))

# # Montagem das Matrizes Globais
# for i in range(Nel):
#     K_temp[4*i:4*i+8, 4*i:4*i+8] = Ke
#     M_temp[4*i:4*i+8, 4*i:4*i+8] = Me
#     G_temp[4*i:4*i+8, 4*i:4*i+8] = Ge
#     Kg += K_temp
#     Mg += M_temp
#     Gg += G_temp

# alpha = 1E-4
# Ce = alpha*Kg

# '''
# Adicionar efeitos do disco no sistema (nó 3)
# '''

# # massa do disco
# Md = np.array([[m,0,0,0],
#               [0,m,0,0],
#               [0,0,Id,0],
#               [0,0,0,Id]])

# # efeito giroscópico
# Gd = np.array([[0,0,0,0],
#               [0,0,0,0],
#               [0,0,0,-Ip],
#               [0,0,Ip,0]])

# Mg[8:12,8:12] += Md
# Gg[8:12,8:12] += Gd

# '''
# Procedimento para plotagem do Diagrama de Campbell
# '''
# n = 4
# '''
#   Adicionar efeitos dos Mancais no sistema (nós 1 e 6)
# '''
# mancal_1 = mancal1(W)

# kyy1 = mancal_1[0]
# kyz1 = mancal_1[1]
# kzy1 = mancal_1[2]
# kzz1 = mancal_1[3]

# cyy1 = mancal_1[4]
# cyz1 = mancal_1[5]
# czy1 = mancal_1[6]
# czz1 = mancal_1[7]

# mancal_2 = mancal2(W)

# kyy2 = mancal_2[0]
# kyz2 = mancal_2[1]
# kzy2 = mancal_2[2]
# kzz2 = mancal_2[3]

# cyy2 = mancal_2[4]
# cyz2 = mancal_2[5]
# czy2 = mancal_2[6]
# czz2 = mancal_2[7]    

# # Matriz de rigidez p/ o mancal 1 (nó 1)
# Km1 = np.array([[kyy1,kyz1,0,0],                  
#                 [kzy1,kzz1,0,0],
#                 [0,0,0,0],
#                 [0,0,0,0]])
# # Matriz de rigidez p/ o mancal 2 (nó 6)
# Km2 = np.array([[kyy2,kyz2,0,0],
#                 [kzy2,kzz2,0,0],
#                 [0,0,0,0],
#                 [0,0,0,0]])

# # Matriz de amortecimento p/ o mancal 1 (nó 1)
# Cm1 = np.array([[cyy1,cyz1,0,0],
#                 [czy1,czz1,0,0],
#                 [0,0,0,0],
#                 [0,0,0,0]])
# # Matriz de amortecimento p/ o mancal 2 (nó 6)
# Cm2 = np.array([[cyy2,cyz2,0,0],
#                 [czy2,czz2,0,0],
#                 [0,0,0,0],
#                 [0,0,0,0]])

# # atualização dos mancais na Matriz de Rigidez Global
# Kg[0:n,0:n] += Km1
# Kg[20:20+n,20:20+n] += Km2

# # atualização dos mancais na Matriz de Amortecimento Global
# Ce[0:n,0:n] += Cm1
# Ce[20:20+n,20:20+n] += Cm2    

# Cg = Ce + W*Gg

# sys = vbt.VibeSystem(Mg,Cg,Kg)


v1 = [9.284108452416268e-10-3.1034089507569e-09j,
      0.0005582870659737103-2.9356267483881807e-06j,
      0.0008864415827315538-4.619546625380575e-06j,
      0.0008410839015670674-4.346664881165814e-06j,
      0.0004995108203054564-2.5750185519104172e-06j,
      8.435433793234114e-10-2.5558091190952847e-09j]

w1 =[6.532831584560326e-10-8.448518377833097e-10j,
     2.9676540030822034e-06+0.0005962199150490342j,
     4.67819529960179e-06+0.000948194369081047j,
     4.407957652252124e-06+0.0009010456360659645j,
     2.612032405728007e-06+0.000535438501766759j,
     7.021828660350036e-10-7.591993029349955e-10j]

v2 = [-8.37923853886091e-06+3.1469839836545475e-05j,
      -0.037937816311087204-0.0004722129816110074j,
      -0.011483372575864858-0.0007806268614068399j,
      0.054307694499350495+0.0012027295966975243j,
      0.07330144381758169+0.005481272563331065j,
      0.03920683037614693+0.00967932269324905j]

w2 = [3.0501591611455034e-05+9.102651053451432e-06j,
      -0.0006126653473616708+0.037135455562401956j,
      -0.0008536251441193703+0.011121288251447054j,
      0.0014984322363470921-0.05306854202440596j,
      0.006101397946081812-0.07144464922891708j,
      0.010340911087601778-0.03849134469081091j]

v22 = [-7.322008909334934e-07-3.855116086394669e-06j,
      0.03480611502871213+5.5636900127618676e-05j,
      0.00804193643024194+1.492031573900772e-05j,
      -0.04813423448443293-0.00037186710475450423j,
      -0.04986346686209195-0.0005442513642372922j,
      9.40927098063127e-07+4.002377011515535e-06j]

w22 = [-3.7775932173156946e-06+5.973354632354435e-07j,
       2.6081940614538363e-05-0.03402041642139457j,
       6.139613181200751e-06-0.0077412259763122455j,
       -0.0003173816189697351+0.04689833060058275j,
       -0.0004771278427055242+0.04832007061805051j,
       3.960074247433286e-06-8.43287486961485e-07j]


v3 = [2.331057234408343e-08-5.036725685308983e-07j,
      0.02880338311922634+0.0014268277745047302j,
      -0.005350097545988522-1.4759552905753014e-05j,
      0.00924001843144997-0.0008218350651430487j,
      0.04510801735038957-0.00023982727155835715j,
      1.346988949563865e-09-7.242334338334366e-07j]

w3 = [-4.724419518675027e-07-4.168495913339169e-08j,
      0.001292016579990641-0.027243519189459164j,
      3.714171416139677e-06+0.005224134794728388j,
      -0.0008292031826936917-0.009860917022148157j,
      -0.00037471118688539324-0.04422892134637933j,
      -7.008702926297842e-07-2.108797436566102e-08j]

X0 = [0,0,0,0,0,0]

X = np.linspace(1,6,NumNodes) # Nodal coordinates
fig1, ax1 = plt.subplots()
ax1.set_xlabel('autovetor')
ax1.set_ylabel('nó')
ax1.plot(X,w3,'b',label='vm')
ax1.plot(X,X0,'orange',label='vm')
plt.savefig('modo3_w.png',dpi=600,bbox_inches = 'tight')