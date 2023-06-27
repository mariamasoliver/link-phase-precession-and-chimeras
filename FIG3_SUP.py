#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:46:14 2023

@author: maria
"""

import matplotlib.pyplot as plt
import numpy as np 
import os


path_files = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files'
os.chdir(path_files)
w = 2.8
data=np.loadtxt('theta_A_0.95_N_500_beta_0.2_ic_chimera_om_2.8_tfinal_10000.txt')
dtfile = 0.1

theta_ring = data[:,1:]
time_ring = data[:,0]


N = len(theta_ring[0,:])
costheta =  np.cos(theta_ring)
lfp_ring = np.mean(costheta, axis=1)

#%%
###Two pop
data1=np.loadtxt("phi_A_0.1_N_3_beta_0.025_ic_0_omega_2.6.txt")
data2=np.loadtxt("theta_A_0.1_N_3_beta_0.025_ic_0_omega_2.6.txt")

theta_twopop = np.append(data1[:,1:],data2[:,1:],axis=1)
time_twopop = data1[:,0]
costheta =  np.cos(theta_twopop)
lfp_twopop = np.mean(costheta, axis=1)
#%%
from datetime import date
lighteal = '#74a9cf'  
darkteal = '#034e7b'
myred = '#bb566b'
params = {
    'axes.labelsize': 28,
        'xtick.labelsize': 28,
            'ytick.labelsize': 28,
                'legend.fontsize': 28,
                'axes.titlesize':28,
                    'text.usetex': False,
                        #'font': 'Helvetica',
                        'mathtext.bf': 'helvetica:bold',

                    }

plt.rcParams.update(params)
fig = plt.figure(figsize = (12,3.25))

ax = fig.add_subplot(121)
ax.plot(time_ring,lfp_ring+2,color = myred, lw = 2)
ax.plot(time_ring,np.cos(theta_ring[:,100])-0.5,color='#ebbd30', lw = 2)
ax.set_xlim(1000,1050)
ax.set_ylim(-1.5,3)
ax.set_yticklabels([])
# ax.axis("off")
ax1 = fig.add_subplot(122)
ax1.plot(time_twopop,lfp_twopop+2,color = darkteal, lw = 2)
ax1.plot(time_twopop,np.cos(theta_twopop[:,3])-0.5,color='#ebbd30', lw = 2)
ax1.set_xlim(1000,1050)
ax1.set_ylim(-1.5,3)
ax1.set_yticklabels([])
# ax1.axis("off")
path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)

plt.tight_layout()
today = date.today()
figure_nom = 'FIG3_sup_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()
