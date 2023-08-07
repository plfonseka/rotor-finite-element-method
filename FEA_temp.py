# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 00:00:17 2020

@author: Lucas
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from FEA_function import FEA

x0 = np.zeros((24,1))
xp0 = np.zeros((24,1))
aux = np.block([[x0],[xp0]])
z0 = np.reshape(aux,(len(aux),))

tin = 0
tf = 2

sol = solve_ivp(FEA,[tin,tf],z0,method='Radau')

T = sol.t
Z = sol.y

dof = len(x0)
x = Z[0:dof,:]
xp = Z[dof:2*dof,:]

i = 0

# ########### PLOT Uy #############
# fig1, ax1 = plt.subplots()
# ax1.set_xlabel('tempo [s]')
# ax1.set_ylabel('Deslocamento [m]')
# ax1.plot(T,x[0,i:],'b',label='Ym')