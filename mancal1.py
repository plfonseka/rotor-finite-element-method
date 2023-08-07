# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 22:32:47 2020

@author: Lucas
"""
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def mancal1():

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
    F1 = 0.6*P
    
    kyy1 = []
    kyz1 = []
    kzy1 = []
    kzz1 = []
    
    cyy1 = []
    cyz1 = []
    czy1 = []
    czz1 = []
    
    W1 = []
    
    for omega in np.arange(1,2400,5):
        
        Fneta = (neta *(lm**3)*omega*rm) / (2*cr**2)
        x0 = 0.8
        
        def S1(eps):
            return ((np.pi/2) * (eps/((1-eps**2)**2)) * (np.sqrt(1-eps**2) + (4*eps/np.pi))) - F1/Fneta
        
        f_eps = fsolve(S1,x0)
        feps = f_eps[0]
    
        Ae = 4 / (np.pi**2 + ((16 - np.pi**2)*feps**2))**(3/2)
    
        aux1 = F1/cr
        
        g11 = (2*np.pi**2 + (16 - np.pi**2) * feps**2) * Ae * aux1
        g12 = (np.pi/4)*((np.pi**2 - 2*np.pi**2 * feps**2 - (16 - np.pi**2)*feps**4 ) / feps*(1 - feps**2)**(1/2)) * Ae * aux1
        g21 = (-np.pi/4)*((np.pi**2 + (32 + np.pi**2)*feps**2 + (32 - 2*np.pi**2)*feps**4)/(feps*(1 - feps**2)**(1/2))) * Ae * aux1
        g22 = ((np.pi**2 + (32 + np.pi**2)*feps**2 + (32 - 2*np.pi**2)*feps**4)/((1 - feps**2)**(1/2))) * Ae * aux1
    
        kyy1.append(g11)
        kyz1.append(g12)
        kzy1.append(g21)
        kzz1.append(g22)
            
        aux2 = F1/(cr*omega)
        
        bt11 = (np.pi/2)*(np.sqrt(1 - feps**2) / feps)*(np.pi**2 + (2*np.pi**2 - 16)*feps**2 ) * Ae * aux2
        bt12 = -(2 * np.pi**2 + (4*np.pi**2 - 32)*feps**2 ) * Ae * aux2
        bt21 = bt12
        bt22 = (np.pi/2)*((np.pi**2 + (48 - 2*np.pi**2)*feps**2 + (np.pi**2 * feps**4)) / (np.sqrt(1 - feps**2) * feps)) * Ae * aux2
            
        cyy1.append(bt11)
        cyz1.append(bt12)
        czy1.append(bt21)
        czz1.append(bt22)
        
        W1.append(omega)
        
    return kyy1,kyz1,kzy1,kzz1,cyy1,cyz1,czy1,czz1