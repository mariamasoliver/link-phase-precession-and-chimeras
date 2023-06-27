#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 16:57:09 2023

@author: maria
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from scipy.signal import find_peaks

path = '/Users/maria/Documents/TempCSG/LAING_setup'
os.chdir(path)



w_vec = np.array([1, 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.5 ,
        1.55, 1.6 , 1.65, 1.7 , 1.75, 1.8 , 1.85, 1.9 , 1.95, 2, 2.05,
        2.1 , 2.15, 2.2 , 2.25, 2.3 , 2.35, 2.4 , 2.45, 2.5 , 2.55, 2.6 ,
        2.65, 2.7 , 2.75, 2.8 , 2.85, 2.9 , 2.95, 3])


# w_vec = np.array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
mpv_un = np.zeros((len(w_vec),2))
mpv_unsyn = np.zeros((len(w_vec),2))
mpv_syn = np.zeros(len(w_vec))


    
n = 0
for omega in w_vec:
    mpv = np.loadtxt('mpv_omega_'+str(omega)+'_N_500.txt')

    sincro = np.min(mpv)
    mask_1 = (mpv == sincro)
    mask_2 = (mpv != sincro)
    mpv_unsyn_vec = mpv[mask_2]
    mpv_syn[n] = sincro
    
    mpv_unsyn[n,0] = np.mean(mpv_unsyn_vec)
    mpv_unsyn[n,1] = np.std(mpv_unsyn_vec) 
    
    n = n +1


N = 10
n = 0
for omega in w_vec:
    data1 = np.loadtxt('uncoupled_theta_A_0.95_N_10_beta_0.2_ic_chimera_om_'+str(omega)+'_tfinal_1000.txt')
    theta = data1[:,1:]
    mpv = np.zeros(N)
    for j in range(N):
        y = theta[:,j]
        peaks_t = find_peaks(y)
        peaks_1 = peaks_t[0]
        rotations_theta = len(peaks_1)
        peaks_1 = peaks_t[0]
        mean_phase_velocity_theta = 2*np.pi*(rotations_theta)/1000
        mpv[j] = np.round(np.mean(mean_phase_velocity_theta),3)
    
    mpv_un[n,0] = np.mean(mpv)
    mpv_un[n,1] = np.std(mpv) 
    n = n+1
    

# for i in range(6):
#     mpv_unsyn[4+i,0] = mpv_syn[4+i]
#     mpv_unsyn[4+i,1] = 0

ti = 1000
tf = 3000
dtfile = 0.1
nti = int(ti/dtfile)
ntf = int(tf/dtfile)


data_u=np.loadtxt('uncoupled_theta_A_0.95_N_10_beta_0.2_ic_chimera_om_2.8_tfinal_1000.txt')
theta_u = data_u[:,1:]
data_u=np.loadtxt('uncoupled_theta_A_0.95_N_10_beta_0.2_ic_chimera_om_2.6_tfinal_1000.txt')
theta_u_bis = data_u[:,1:]
data=np.loadtxt('theta_A_0.95_N_500_beta_0.2_ic_chimera_om_2.8_tfinal_10000.txt',skiprows = nti, max_rows = ntf)
dtfile = 0.1
#we use the same delta we used to compute the mpv
deltat = 1000
ndeltat = int(deltat/dtfile)
theta_ring = data[ndeltat:2*ndeltat,1:]
time_ring = data[ndeltat:2*ndeltat,0]



ring_mpv = np.loadtxt('mpv_omega_2.8_N_500.txt') 

nodes = 10
mpv_max = np.max(ring_mpv)
mpv_min = np.min(ring_mpv)

mpv_max_node = np.where(ring_mpv==mpv_max)
mpv_min_node = np.where(ring_mpv==mpv_min)

unsin = theta_ring[:,mpv_max_node[0][0]:mpv_max_node[0][0]+nodes]
sin = theta_ring[:,mpv_min_node[0]]


#%%
path = '/Users/maria/Documents/TempCSG/twopop_setup'
os.chdir(path)

mpv_unsyn_tp = np.zeros(len(w_vec))
mpv_syn_tp = np.zeros(len(w_vec))

n = 0
for omega in w_vec:
    mpv = np.loadtxt('mpv_A_0.1_delta_1000.0_ic_0_w_'+str(omega)+'_TC_1_N_3.txt')

    sincro = np.min(mpv)
    mask_1 = (mpv == sincro)
    mask_2 = (mpv != sincro)
    mpv_unsyn_tp[n]= mpv[mask_2]
    mpv_syn_tp[n] = sincro
    
    n = n +1
    
path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files/two_pop_chimera'
os.chdir(path)

w = 2.6
data1=np.loadtxt('theta_A_0.1_N_3_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)
data2=np.loadtxt('phi_A_0.1_N_3_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)

theta_tp = data1[ndeltat:2*ndeltat,1:]
phi_tp = data2[ndeltat:2*ndeltat,1:]

time_tp = data1[ndeltat:2*ndeltat,0]

#%%

params = {
    'axes.labelsize': 22,
        'xtick.labelsize': 22,
            'ytick.labelsize': 22,
                'legend.fontsize': 22,
                'axes.titlesize':22,
                    'text.usetex': False,
                        #'font': 'Helvetica',
                        'mathtext.bf': 'helvetica:bold',

                    }

plt.rcParams.update(params)
fig = plt.figure(figsize = (20,7))

ax1 = fig.add_subplot(2,2,1)

lighteal = '#74a9cf'  
darkteal = '#034e7b'
myred = '#bb566b'

ax1.plot(w_vec, mpv_unsyn[:,0], marker = 'o', lw = '1.5', markersize=7, color = myred, alpha=0.6)
plt.fill_between(w_vec, mpv_unsyn[:,0] - mpv_unsyn[:,1], mpv_unsyn[:,0] + mpv_unsyn[:,1],color=myred, alpha=0.4)
ax1.plot(w_vec, mpv_syn, marker = '^', lw = '2', markersize=7, color = myred)
ax1.plot(w_vec, mpv_un[:,0], marker = 's', lw = '2', markersize=7, color = 'grey')
ax1.set_ylim(0,3.2)

plt.xlabel(r'Intrinsic frequency $\rho$')
plt.ylabel('Mean phase\n velocity $\Omega$')
# ax1.set_ylim(0.25,1.1)

ax3 = fig.add_subplot(2,2,3)
   
ax3.plot(time_ring, unsin[:,3], color=myred,alpha = 0.6, linewidth=4)
ax3.plot(time_ring, sin[:,0]-8, color=myred, linewidth=4)
ax3.plot(time_ring, theta_u[:,0]+8, color='grey', linewidth=4)

ax3.set_xlim(2100,2140)
ax3.set_xlabel('time $t$') 


ax2 = fig.add_subplot(2,2,2)

lighteal = '#74a9cf'  
darkteal = '#034e7b'
myred = '#bb566b'

ax2.plot(w_vec, mpv_unsyn_tp, marker = 'o', lw = '1.5', markersize=7, color = darkteal, alpha=0.6)
ax2.plot(w_vec, mpv_syn_tp, marker = '^', lw = '2', markersize=7, color = darkteal)
ax2.plot(w_vec, mpv_un[:,0], marker = 's', lw = '2', markersize=7, color = 'grey')
ax2.set_ylim(0,3.2)

plt.xlabel(r'Intrinsic frequency $\rho$')
# plt.ylabel('Mean phase\n velocity $\Omega$')
# ax1.set_ylim(0.25,1.1)

ax3 = fig.add_subplot(2,2,4)
   
ax3.plot(time_tp, phi_tp[:,0], color=darkteal,alpha = 0.6, linewidth=4)
ax3.plot(time_tp, theta_tp[:,0]-8, color=darkteal, linewidth=4)
ax3.plot(time_tp, theta_u_bis[:,0]+8, color='grey', linewidth=4)

ax3.set_xlim(2100,2140)
ax3.set_xlabel('time $t$') 

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)

plt.tight_layout()
today = date.today()
figure_nom = 'FIG6_sup_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()
