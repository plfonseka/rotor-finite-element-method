# -*- coding: utf-8 -*-
"""
Universidade Estadual de Campinas
IM 342 - Análise de Máquinas Rotativas
Pedro Lucas - Ra 263117

Resposta em Frequência do Sistema
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt  
# import vibration_toolbox as vbt

def plot_freq_resp(self, modes=None, ax0=None, ax1=None, **kwargs):
    
    if ax0 is None or ax1 is None:
        fig, ax = plt.subplots(2)
        if ax0 is not None:
            _, ax1 = ax
        if ax1 is not None:
            ax0, _ = ax
        else:
            ax0, ax1 = ax

    omega, magdb, phase = self.freq_response(modes=modes)


    # ax0.plot(omega, magdb[20,20], **kwargs)
    # ax1.plot(omega, phase[20,20], **kwargs)

    # ax0.plot(omega, magdb[20,21], **kwargs)
    # ax1.plot(omega, phase[20,21], **kwargs)
    
    # ax0.plot(omega, magdb[21,21], **kwargs)
    # ax1.plot(omega, phase[21,21], **kwargs)

    # ax0.plot(omega, magdb[21,20], **kwargs)
    # ax1.plot(omega, phase[21,20], **kwargs)
    
    # for i in range(20,22):        
    #     ax0.plot(omega, magdb[0,i], **kwargs)
    #     ax1.plot(omega, phase[0,i], **kwargs)

    # for i in range(0,2):
    #     ax0.plot(omega, magdb[0,i], **kwargs)
    #     ax1.plot(omega, phase[0,i], **kwargs)  
    
    for i in range(8,10):
        ax0.plot(omega, magdb[0,i], **kwargs)
        ax1.plot(omega, phase[0,i], **kwargs)  
    
    for ax in [ax0, ax1]:
        ax.set_xlim(0,max(omega))
        ax.yaxis.set_major_locator(
            mpl.ticker.MaxNLocator(prune='lower'))
        ax.yaxis.set_major_locator(
            mpl.ticker.MaxNLocator(prune='upper'))
    
    ax0.set_title('Resposta em Frequência do Sistema Amortecido')
    ax0.set_ylabel('Magnitude $(dB)$')
    # legend = ax0.legend(loc='upper right', shadow = True)
    # legend.get_frame().set_facecolor('white')
    ax1.set_ylabel('Ângulo de Fase $(°)$')
    ax1.set_xlabel('Frequência (rad/s)')
        
    return ax0, ax1