#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  17 14:30:35 2023

@author: maria
"""

import matplotlib.pyplot as plt
import numpy as np 
import os
import matplotlib.cm as cm
from datetime import date

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files/two_pop_chimera'
os.chdir(path)


"""
Loading data 
RATIO MPV
"""
meanratio1 = np.loadtxt('ratio_mpv_A_0.1_delta_1000.0_ic_0_TC_1_N_25_two_pop.txt')

ti = 1000
tf = 9000
dtfile = 0.1
nti = int(ti/dtfile)
ntf = int(tf/dtfile)
nodes = 10

"""
w = 1.6
"""
w = 1.6
data1=np.loadtxt('theta_A_0.1_N_25_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)
data2=np.loadtxt('phi_A_0.1_N_25_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)

deltat = 1000
ndeltat = int(deltat/dtfile)
#we use the same delta we used to compute the mpv
theta_w1 = data1[:,1:]
phi_w1 = data2[:,1:]

time = np.arange(0,len(theta_w1),dtfile)

"""1. Poincare map"""

unsin = theta_w1
sin = phi_w1


spikes_un_w1 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(unsin)-1):
        if (unsin[i,j] < 0 and unsin[i+1,j] > 0):
            spikes = np.append(spikes,time[i])
    spikes_un_w1.append(spikes)

spikes_sin_w1 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(sin)-1):
        if (sin[i,j] < 0 and sin[i+1,j] > 0):
            spikes = np.append(spikes,time[i])
    spikes_sin_w1.append(spikes)
    
    
"""
Rotations
"""
from scipy.signal import find_peaks

dtplot = 0.1
deltat = 1000
ntf_output = len(phi_w1[:,0])
nodes_tot = len(phi_w1[0,:])*2
ndeltat = int(deltat/dtplot)
ntotal_deltat = int(ntf_output/ndeltat)
ndeltat_vec = np.arange(ntotal_deltat)*ndeltat

phase = np.append(theta_w1,phi_w1,axis=1)

j = 0
rotations_phase = np.zeros((ntotal_deltat,nodes_tot))
for l in range (nodes_tot):
    peaks_t = find_peaks(phase[:,l], height=2.5)
    peaks_1 = peaks_t[0]
    for i in range (len(ndeltat_vec)):
        rotations_phase[i, l] = len(np.nonzero((peaks_1<ndeltat_vec[i])*1)[0]) - sum(rotations_phase[:i, l])

      
"""
MPV
"""

mpv = np.zeros((ntotal_deltat-1,nodes_tot))
for i in range(nodes_tot):
    mpv[:,i] = (2*np.pi*(rotations_phase[1:,i])/deltat)
    
mean_mpv_w1 = np.mean(mpv, axis=0)


mask_s = (mean_mpv_w1==min(mean_mpv_w1))
mean_mpv_w1_sin = mean_mpv_w1[mask_s]  

mean_mpv_w1_sin = np.round(mean_mpv_w1_sin,2)  
  
mask_us = (mean_mpv_w1 != min(mean_mpv_w1))
mean_mpv_w1_usin = mean_mpv_w1[mask_us]
    
mean_mpv_w1_usin = np.round(mean_mpv_w1_usin,2)  

"""
w = 2.6
"""
w = 2.6
data1=np.loadtxt('theta_A_0.1_N_25_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)
data2=np.loadtxt('phi_A_0.1_N_25_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)

theta_w2 = data1[:,1:]
phi_w2 = data2[:,1:]

time = np.arange(0,len(theta_w1),dtfile)

"""1. Poincare map"""

unsin = theta_w2
sin = phi_w2

spikes_un_w2 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(unsin)-1):
        if (unsin[i,j] < 0 and unsin[i+1,j] > 0):
            spikes = np.append(spikes,time[i])
    spikes_un_w2.append(spikes)

spikes_sin_w2 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(sin)-1):
        if (sin[i,j] < 0 and sin[i+1,j] > 0):
            spikes = np.append(spikes,time[i])
    spikes_sin_w2.append(spikes)
    
dtplot = 0.1
deltat = 1000
ntf_output = len(phi_w2[:,0])
nodes_tot = len(phi_w2[0,:])*2
ndeltat = int(deltat/dtplot)
ntotal_deltat = int(ntf_output/ndeltat)
ndeltat_vec = np.arange(ntotal_deltat)*ndeltat

phase = np.append(theta_w2,phi_w2,axis=1)

j = 0
rotations_phase = np.zeros((ntotal_deltat,nodes_tot))
for l in range (nodes_tot):
    peaks_t = find_peaks(phase[:,l], height=2.5)
    peaks_1 = peaks_t[0]
    for i in range (len(ndeltat_vec)):
        rotations_phase[i, l] = len(np.nonzero((peaks_1<ndeltat_vec[i])*1)[0]) - sum(rotations_phase[:i, l])

     
"""
MPV
"""

mpv = np.zeros((ntotal_deltat-1,nodes_tot))
for i in range(nodes_tot):
    mpv[:,i] = (2*np.pi*(rotations_phase[1:,i])/deltat)
    
mean_mpv_w2 = np.mean(mpv, axis=0)

mask_s = (mean_mpv_w2==min(mean_mpv_w2))
mean_mpv_w2_sin = mean_mpv_w2[mask_s]  
mean_mpv_w2_sin = np.round(mean_mpv_w2_sin,2)
    
mask_us = (mean_mpv_w2 != min(mean_mpv_w2))
mean_mpv_w2_usin = mean_mpv_w2[mask_us] 
mean_mpv_w2_usin = np.round(mean_mpv_w2_usin,2)
#%%

from datetime import date

"""
plot Ratio sin/unsin 1pop
"""


w_vec = np.array([1, 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.5 ,
       1.55, 1.6 , 1.65, 1.7 , 1.75, 1.8 , 1.85, 1.9 , 1.95, 2, 2.05,
       2.1 , 2.15, 2.2 , 2.25, 2.3 , 2.35, 2.4 , 2.45, 2.5 , 2.55, 2.6 ,
       2.65, 2.7 , 2.75, 2.8 , 2.85, 2.9 , 2.95, 3])

params = {
    'axes.labelsize': 18,
        'xtick.labelsize': 18,
            'ytick.labelsize': 18,
                'legend.fontsize': 18,
                    
                    'text.usetex': False,
                        #'font': 'Helvetica',
                        'mathtext.bf': 'helvetica:bold',

                    }

plt.rcParams.update(params)

fig = plt.figure(figsize = (7.5,10))
ax1 = fig.add_subplot(3,2,(5,6))



myred = '#bb566b'
mylightred = '#D699A6'

plt.plot(w_vec, meanratio1, marker = 'o', markersize=6, color = 'k')
#plt.fill_between(w_vec, meanratio1[:,1] - meanratio1[:,2], meanratio1[:,1] + meanratio1[:,2],color=myred, alpha=0.4)
ax1.set_ylabel(r'$\Omega_{ratio}$', size = '18')
plt.xlabel(r'$\rho$', size = '18')
plt.yticks()
plt.xticks()

# plt.axvline(2.8,lw = 10,alpha=0.4, color='grey')
# plt.axhline(0.88,lw = 10,alpha=0.4, color='grey')
plt.xlabel(r'$\rho$', size = '18')
ax1.set_ylim(0,1)


"""
poncaire maps
general info
"""

num = nodes
colors = []
cmap = cm.get_cmap('Blues')
step = (1-0)/25
for i in range(25):
    ind = 0.3+i*step
    colors.append(cmap(ind))
    
myred = '#bb566b'
mylightred = '#D699A6'
mygrey = '#969696'

darkgreen = '#225522'
lightgreen = '#90AA90'

"""
poncaire map
w = 1.8
"""
ax2 = fig.add_subplot(3,2,1)

nspikes = 30
for s in range(2,nodes):
    x_un = spikes_un_w1[s]
    x_sin = spikes_sin_w1[s]
    ax2.plot(x_un[:nspikes], np.zeros(nspikes)+0.5*(s+3)+0.25, '|', markersize=18, markeredgewidth=4, markeredgecolor=colors[s])
    # ax2.plot(x_sin[:nspikes], np.zeros(nspikes)-3.75, '|', markersize=18, markeredgewidth=4, markeredgecolor=mygrey)
# for p in range(nspikes):
    # plt.vlines(x_sin[p], ymin=-4.5,ymax= 0, linestyle = 'dotted', color = 'grey',linewidth=2)
# ax2.set_ylim(-4,0)
ax2.set_xlim(0,30)

# ax2.axis("off")

"""
mpv
w = 1.8
"""
node_vec = np.arange(0,len(mean_mpv_w1_usin))
ax2b = fig.add_subplot(3,2,3)

ax2b.plot(np.arange(len(mean_mpv_w1_sin))+len(mean_mpv_w1_sin), mean_mpv_w1_sin, marker = 'o', markersize=4, linewidth=0, color = 'k')
for n in range (len(mean_mpv_w2_sin)):
    ax2b.plot(node_vec[n], mean_mpv_w1_usin[n], marker = 'o', markersize=4, linewidth=0, color = 'k')

ax2b.set_ylabel(r'$\Omega$', size = '18')
plt.xlabel(r'Oscillator $i$', size = '18')
# ax2b.set_ylim(-2.5,0)
# ax2b.set_xlim(0,25)
# ax2b.axis("off")

"""
poncaire map
w = 2.8
"""
ax3 = fig.add_subplot(3,2,2)


nspikes = 30  
for s in range(3,nodes):
    x_un = spikes_un_w2[s]
    x_sin = spikes_sin_w2[s]
    ax3.plot(x_un[:nspikes], np.zeros(nspikes)-0.5*(s-3)-0.25, '|', markersize=18, markeredgewidth=4, markeredgecolor=colors[s])
    ax3.plot(x_sin[:nspikes], np.zeros(nspikes)-3.75, '|', markersize=18, markeredgewidth=4, markeredgecolor=mygrey)
for p in range(nspikes):
    plt.vlines(x_sin[p], ymin=-4.5,ymax= 0, linestyle = 'dotted', color = 'grey',linewidth=2)
ax3.set_ylim(-4,0)
ax3.set_xlim(0,30)
# ax3.axis("off")

"""
mpv
w = 2.8
"""
node_vec = np.arange(0,len(mean_mpv_w2_usin))
ax3b = fig.add_subplot(3,2,4)
ax3b.plot(np.arange(len(mean_mpv_w2_sin))+len(mean_mpv_w2_sin), mean_mpv_w2_sin, marker = 'o', markersize=4, linewidth=0, color = 'k')
for n in range (len(mean_mpv_w2_sin)):
    ax3b.plot(node_vec[n], mean_mpv_w2_usin[n], marker = 'o', markersize=4, linewidth=0, color = 'k')

plt.xlabel(r'Oscillator $i$', size = '18')
# ax2b.set_ylim(-2.5,0)
# ax2b.set_xlim(0,25)
# ax2b.axis("off")


plt.tight_layout()
today = date.today()

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)

figure_nom = 'FIG5_suplementary_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()