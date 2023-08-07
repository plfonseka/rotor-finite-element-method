# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 16:18:46 2020

@author: Lucas
"""
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D



t1 = open("t3.txt").read().split()
v1 = open("mancal1_80vm.txt").read().split()
w1 = open("mancal1_80wm.txt").read().split()

t2 = open("t3.txt").read().split()
v2 = open("disco_80vm.txt").read().split()
w2 = open("disco_80wm.txt").read().split()

t3 = open("t3.txt").read().split()
v3 = open("mancal2_80vm.txt").read().split()
w3 = open("mancal2_80wm.txt").read().split()

a1 = np.loadtxt(t1,dtype='float')
b1 = np.loadtxt(v1,dtype='float')
c1 = np.loadtxt(w1,dtype='float')

a2 = np.loadtxt(t2,dtype='float')
b2 = np.loadtxt(v2,dtype='float')
c2 = np.loadtxt(w2,dtype='float')

a3 = np.loadtxt(t3,dtype='float')
b3 = np.loadtxt(v3,dtype='float')
c3 = np.loadtxt(w3,dtype='float')

# X = np.linspace(0,5,7405)


# fig = plt.figure()
# ax = plt.axes(projection='3d')

# z_line = c1
# x_line = X
# y_line = b1
# ax.plot3D(x_line, y_line, z_line, 'k')


######### PLOT mancal1 #############
fig1, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(a1,b1,'b',label='vm')
plt.savefig('m1_80vm.png',dpi=600,bbox_inches = 'tight')

fig2, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(a1,c1,'b',label='wm')
plt.savefig('m1_80wm.png',dpi=600,bbox_inches = 'tight')

########## PLOT disco #############
fig3, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(a2,b2,'b',label='wm')
plt.savefig('disco_80vm.png',dpi=600,bbox_inches = 'tight')

fig4, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(a2,c2,'b',label='wm')
plt.savefig('disco_80wm.png',dpi=600,bbox_inches = 'tight')

########## PLOT mancal2 #############
fig5, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(a3,b3,'b',label='wm')
plt.savefig('m2_80vm.png',dpi=600,bbox_inches = 'tight')

fig6, ax1 = plt.subplots()
ax1.set_xlabel('tempo [s]')
ax1.set_ylabel('Deslocamento [m]')
ax1.plot(a3,c3,'b',label='wm')
plt.savefig('m2_80wm.png',dpi=600,bbox_inches = 'tight')
