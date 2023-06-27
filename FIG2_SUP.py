#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:38:11 2023

@author: maria
"""


import matplotlib.pyplot as plt
import numpy as np 
import os
from datetime import date

path_files = '/Users/maria/Documents/TempCSG/with_phase_perturbation_vs_2/'
os.chdir(path_files)
epsilon = 0.5
perturbation=np.loadtxt('random_kick_time_A_0.1_N_3_beta_0.025_seed_4_epsilon_'+str(epsilon)+'_omega_2.6_skip_nt_500000_1_kick.txt')
theta = np.loadtxt('theta_A_0.1_N_3_beta_0.025_seed_4_epsilon_'+str(epsilon)+'_omega_2.6_skip_nt_500000_1_kick.txt')
phi=np.loadtxt('phi_A_0.1_N_3_beta_0.025_seed_4_epsilon_'+str(epsilon)+'_omega_2.6_skip_nt_500000_1_kick.txt')
dtfile = 0.1

time = theta[:,0]
theta = theta[:,1:]
phi = phi[:,1:]

#%%
def R_v1(spacetime):
    """Returns the Kuramoto order parameter for an TxN numpy array.
    Parameters:
    -----------
    spacetime: TxN numpy array (N: number of nodes, T: time).
               Each row represents a timeseries of the corresponding node.
    Returns:
    --------
    r: Global order parameter (Scalar number at each time)
    """
    scos = np.cos(spacetime).sum(axis=1)
    ssin = np.sin(spacetime).sum(axis=1)
    r = np.sqrt((scos*scos + ssin*ssin)) / (1.0*spacetime.shape[1])
    return r  

#%%
rvtheta = R_v1(theta)
rvphi = R_v1(phi)

#%%
lighteal = '#74a9cf'    
teal = '#3690c0'
darkteal = '#034e7b'
mygrey = '#969696'

params = {
    'axes.labelsize': 30,
        'xtick.labelsize': 30,
            'ytick.labelsize': 30,
                'legend.fontsize': 30,
                'axes.titlesize':30,
                    'text.usetex': False,
                        #'font': 'Helvetica',
                        'mathtext.bf': 'helvetica:bold',

                    }

plt.rcParams.update(params)
fig = plt.figure(figsize = (20,15))


"""
132 Time-series: 2 pop
"""
ax5 = fig.add_subplot(3,3,(1,3))
ax5.plot(time, rvtheta, color=darkteal, linewidth=4)
ax5.plot(time, rvphi+2, color=teal, linewidth=4)
plt.axvline(perturbation,color='#D41159', alpha = 0.8, linestyle="--",linewidth=4)
plt.axvline(perturbation+700-10,color='k', alpha = 1, linestyle="-",linewidth=2)
plt.axvline(perturbation-10-50,color='k', alpha = 1, linestyle="-",linewidth=2)
plt.axvline(perturbation-10,color='k', alpha = 1, linestyle="-",linewidth=2)
plt.axvline(perturbation+700-50-10,color='k', alpha = 1, linestyle="-",linewidth=2)
ax5.set_xlim(perturbation-100,perturbation+700)
ax5.plot(perturbation,1,'or')
ax5.set_yticks([0,1,2,3])
ax5.set_yticklabels([0,1,0,1])

ax5.set_ylabel('R(t)')
ax5.set_xlabel('Time (arb. unit)')

ax2 = fig.add_subplot(3,3,4)
ax2.plot(time, rvtheta, color=darkteal, linewidth=4)
ax2.plot(time, rvphi+2, color=teal, linewidth=4)
plt.axvline(perturbation,color='#D41159', alpha = 0.8, linestyle="--",linewidth=4)
ax2.set_xlim(perturbation-50-10,perturbation-10)
ax2.plot(perturbation,1,'or')
ax2.set_yticks([0,1,2,3])
ax2.set_yticklabels([0,1,0,1])

ax2.set_ylabel('R(t)')
ax2.set_xlabel('Time (arb. unit)')

ax2 = fig.add_subplot(3,3,6)
ax2.plot(time, rvtheta, color=darkteal, linewidth=4)
ax2.plot(time, rvphi+2, color=teal, linewidth=4)
plt.axvline(perturbation,color='#D41159', alpha = 0.8, linestyle="--",linewidth=4)
# ax5.set_xlim(perturbation[0]-100,perturbation[1]+100)
# ax5.plot(perturbation,np.zeros(len(perturbation))+1,'or')
ax2.set_xlim(perturbation+700-50-10,perturbation+700-10)
ax2.plot(perturbation,1,'or')
ax2.set_yticks([0,1,2,3])
ax2.set_yticklabels([0,1,0,1])

ax2.set_ylabel('R(t)')
ax2.set_xlabel('Time (arb. unit)')
ax = fig.add_subplot(337)
i = 0
for j in range (0,3):  
    print(j)
    ax.plot(time, theta[:,j]+10*(i), color=darkteal, linewidth=4)
    ax.plot(time, phi[:,j]+10*(i+3), color=teal, linewidth=4)
    i = i +1
# ax.set_xlim(perturbation[0]-65,perturbation[0]-10)
ax.set_xlabel('Time (arb. unit)')
ax.set_yticks([-np.pi,np.pi,-np.pi+30,np.pi+30])
ax.set_yticklabels(['0', r'2$\pi$','0', r'2$\pi$'])
ax.set_xlim(perturbation-50-10,perturbation-10)

ax2 = fig.add_subplot(339)
i = 0
for j in range (0,3):  
    print(j)
    ax2.plot(time, theta[:,j]+10*(i), color=darkteal, linewidth=4)
    ax2.plot(time, phi[:,j]+10*(i+3), color=teal, linewidth=4)
    i = i +1
# ax2.set_xlim(perturbation[3]+10,perturbation[3]+65)
ax2.set_xlabel('Time (arb. unit)')
ax2.set_yticks([-np.pi,np.pi,-np.pi+30,np.pi+30])
ax2.set_yticklabels(['0', r'2$\pi$','0', r'2$\pi$'])
ax2.set_xlim(perturbation+700-50-10,perturbation+700-10)

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)
today = date.today()
seed = 4
figure_nom = 'FIG7_SUP'+str(today)+'_epsilon_'+str(epsilon)+'_bis_1_kick.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()

