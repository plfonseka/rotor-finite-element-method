# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 23:17:38 2020

@author: Lucas
"""
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def mancal2(omega):
    gv = 9.8
    dd = 90E-3
    ld = 15E-3
    rho = 7850
    vol_d = (np.pi*((dd/2)**2))*ld
    m = rho*vol_d
    P = m*gv
    
    dm = 30E-3
    rm = dm/2
    lm = 20E-3
    cr = 90E-6
    neta = 0.07
    F2 = 0.4*P
        
    Fneta = (neta *(lm**3)*omega*rm) / (2*cr**2)
    x0 = 0.8
        
    def S2(eps):
        return ((np.pi/2) * (eps/((1-eps**2)**2)) * (np.sqrt(1-eps**2) + (4*eps/np.pi))) - F2/Fneta
    
    f_ep = fsolve(S2,x0)
    fep = f_ep[0]

    Ae = 4 / (np.pi**2 + ((16 - np.pi**2)*fep**2))**(3/2)

    aux3 = F2/cr
    
    h11 = (2*np.pi**2 + (16 - np.pi**2) * fep**2) * Ae * aux3
    h12 = (np.pi/4)*((np.pi**2 - 2*np.pi**2 * fep**2 - (16 - np.pi**2)*fep**4 ) / fep*(1 - fep**2)**(1/2)) * Ae * aux3
    h21 = (-np.pi/4)*((np.pi**2 + (32 + np.pi**2)*fep**2 + (32 - 2*np.pi**2)*fep**4)/(fep*(1 - fep**2)**(1/2))) * Ae * aux3
    h22 = ((np.pi**2 + (32 + np.pi**2)*fep**2 + (32 - 2*np.pi**2)*fep**4)/((1 - fep**2)**(1/2))) * Ae * aux3

    kyy2=h11
    kyz2=h12
    kzy2=h21
    kzz2=h22
        
    aux4 = F2/(cr*omega)
    
    b11 = (np.pi/2)*(np.sqrt(1 - fep**2) / fep)*(np.pi**2 + (2*np.pi**2 - 16)*fep**2 ) * Ae * aux4
    b12 = -(2 * np.pi**2 + (4*np.pi**2 - 32)*fep**2 ) * Ae * aux4
    b21 = b12
    b22 = (np.pi/2)*((np.pi**2 + (48 - 2*np.pi**2)*fep**2 + (np.pi**2 * fep**4)) / (np.sqrt(1 - fep**2) * fep)) * Ae * aux4
        
    cyy2=b11
    cyz2=b12
    czy2=b21
    czz2=b22

    return kyy2,kyz2,kzy2,kzz2,cyy2,cyz2,czy2,czz2