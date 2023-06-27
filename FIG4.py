#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:30:35 2023

@author: maria
"""

import matplotlib.pyplot as plt
import numpy as np 
import os
import matplotlib.cm as cm
from datetime import date

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files'
os.chdir(path)


"""
Loading data 
RATIO MPV
"""
meanratio1 = np.loadtxt('meanratio_vs_w_N_500_different_delta.txt')
#meanratio2 = np.loadtxt('meanratio_vs_w_N_500_one_delta.txt')
#%%
ti = 1000
tf = 3000
dtfile = 0.1
nti = int(ti/dtfile)
ntf = int(tf/dtfile)
nodes = 8

"""
w = 1
"""
w = 1.8
data=np.loadtxt('theta_A_0.95_N_500_beta_0.2_ic_chimera_om_'+str(w)+'_tfinal_10000.txt',skiprows = nti, max_rows = ntf)

deltat = 1000
ndeltat = int(deltat/dtfile)
#we use the same delta we used to compute the mpv
theta_w1 = data[ndeltat:2*ndeltat,1:]
time_ring = np.arange(0,len(theta_w1),dtfile)

"""1. Poincare map"""

#need to know which are the sincro/unsincro neurons
#during the period of time t1 = 0 and t2 = 1000
#we'll use the mpv. The first 10 nodes with the lowest mpv: sincro pop
#the last 10 nodes with the highest mpv: unsincro pop
#the problem with this method is that the nodes will not be ordered regarding 
# their topological position

ring_mpv_w1 = np.loadtxt('mpv_omega_'+str(w)+'_N_500.txt') 
sorted_index_mpv_w1 = np.argsort(ring_mpv_w1)
#unsin = theta_w1[:,sorted_index_mpv_w1[-nodes:]]
unsin = theta_w1[:,-nodes:]
sin = theta_w1[:,sorted_index_mpv_w1[:nodes]]


spikes_un_w1 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(unsin)-1):
        if (unsin[i,j] < 0 and unsin[i+1,j] > 0):
            spikes = np.append(spikes,time_ring[i])
    spikes_un_w1.append(spikes)

spikes_sin_w1 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(sin)-1):
        if (sin[i,j] < 0 and sin[i+1,j] > 0):
            spikes = np.append(spikes,time_ring[i])
    spikes_sin_w1.append(spikes)

"""
w = 2
"""
w = 2.0
data=np.loadtxt('theta_A_0.95_N_500_beta_0.2_ic_chimera_om_'+str(w)+'_tfinal_10000.txt',skiprows = nti, max_rows = ntf)
dtfile = 0.1

#we use the same delta we used to compute the mpv
theta_w2 = data[ndeltat:2*ndeltat,1:]
time_ring = np.arange(0,len(theta_w2),dtfile)

"""1. Poincare map"""

#need to know which are the sincro/unsincro neurons
#during the period of time t1 = 0 and t2 = 1000
#we'll use the mpv. The first 10 nodes with the lowest mpv: sincro pop
#the last 10 nodes with the highest mpv: unsincro pop

ring_mpv_w2 = np.loadtxt('mpv_omega_'+str(w)+'_N_500.txt') 
sorted_index_mpv_w2 = np.argsort(ring_mpv_w2)
# unsin = theta_w2[:,sorted_index_mpv_w2[-nodes:]]
unsin = theta_w2[:,-nodes:]
sin = theta_w2[:,sorted_index_mpv_w2[:nodes]]


spikes_un_w2 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(unsin)-1):
        if (unsin[i,j] < 0 and unsin[i+1,j] > 0):
            spikes = np.append(spikes,time_ring[i])
    spikes_un_w2.append(spikes)

spikes_sin_w2 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(sin)-1):
        if (sin[i,j] < 0 and sin[i+1,j] > 0):
            spikes = np.append(spikes,time_ring[i])
    spikes_sin_w2.append(spikes)
    
"""
w = 2.8
"""
w = 2.8
data=np.loadtxt('theta_A_0.95_N_500_beta_0.2_ic_chimera_om_'+str(w)+'_tfinal_10000.txt',skiprows = nti, max_rows = ntf)
dtfile = 0.1

#we use the same delta we used to compute the mpv
theta_w3 = data[ndeltat:2*ndeltat,1:]
time_ring = np.arange(0,len(theta_w3),dtfile)

"""1. Poincare map"""

#need to know which are the sincro/unsincro neurons
#during the period of time t1 = 0 and t2 = 1000
#we'll use the mpv. The first 10 nodes with the lowest mpv: sincro pop
#the last 10 nodes with the highest mpv: unsincro pop

ring_mpv_w3 = np.loadtxt('mpv_omega_'+str(w)+'_N_500.txt')  
sorted_index_mpv_w3 = np.argsort(ring_mpv_w3)
# unsin = theta_w3[:,sorted_index_mpv_w3[-nodes:]]
unsin = theta_w3[:,:nodes]#I know from the mpv plot that the first 10 nodes are unsyn
sin = theta_w3[:,sorted_index_mpv_w3[:nodes]]


spikes_un_w3 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(unsin)-1):
        if (unsin[i,j] < 0 and unsin[i+1,j] > 0):
            spikes = np.append(spikes,time_ring[i])
    spikes_un_w3.append(spikes)

spikes_sin_w3 = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(sin)-1):
        if (sin[i,j] < 0 and sin[i+1,j] > 0):
            spikes = np.append(spikes,time_ring[i])
    spikes_sin_w3.append(spikes)


#%%

from datetime import date

"""
plot Ratio sin/unsin 1pop
"""


w_vec = np.array([1, 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.5 ,
       1.55, 1.6 , 1.65, 1.7 , 1.75, 1.8 , 1.85, 1.9 , 1.95, 2, 2.05,
       2.1 , 2.15, 2.2 , 2.25, 2.3 , 2.35, 2.4 , 2.45, 2.5 , 2.55, 2.6 ,
       2.65, 2.7 , 2.75, 2.8 , 2.85, 2.9 , 2.95, 3,3.05, 3.1 , 3.15, 3.2 , 3.25, 3.3 , 3.35, 3.4 , 3.45, 3.5])

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

fig = plt.figure(figsize = (10,10))
ax1 = fig.add_subplot(4,2,(7,8))



myred = '#bb566b'
mylightred = '#D699A6'

plt.plot(w_vec, meanratio1[:,1], marker = 'o', markersize=6, color = 'k')
plt.fill_between(w_vec, meanratio1[:,1] - meanratio1[:,2], meanratio1[:,1] + meanratio1[:,2],color=myred, alpha=0.4)
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
cmap = cm.get_cmap('Reds')
step = (0.8-0.2)/6
for i in range(num):
    ind = 0.2+i*step
    colors.append(cmap(ind))
    
myred = '#bb566b'
mylightred = '#D699A6'
mygrey = '#969696'

"""
poncaire map
w = 1.8
"""
ax2 = fig.add_subplot(4,2,(1,3))


nspikes = 25   
for s in range(num):
    x_un = spikes_un_w1[s]
    x_sin = spikes_sin_w1[s]
    ax2.plot(x_un[:nspikes], np.zeros(nspikes)+0.25*s+0.15, '|', markersize=18, markeredgewidth=4, markeredgecolor=colors[s])
for s in range(3):
    ax2.plot(x_sin[:nspikes], np.zeros(nspikes)-0.25*s-0.5, '|', markersize=18, markeredgewidth=4, markeredgecolor=mygrey)

for p in range(nspikes):
    plt.axvline(x_sin[p],linestyle = 'dotted', color = 'grey',linewidth=2)
# s = 0
# x_un = spikes_un_w1[s]
# ax2.plot(x_un[:nspikes], np.zeros(nspikes)-0.25*s-0.15, '|', markersize=18, markeredgewidth=4, markeredgecolor=myred)
# ax2.set_ylim(-2.5,0)
ax2.set_xlim(0,25)
# ax2.axis("off")

"""
mpv
w = 1.8
"""
ax2b = fig.add_subplot(4,2,5)

plt.plot(np.arange(len(ring_mpv_w1)), ring_mpv_w1, marker = 'o', markersize=3, color = 'k')
ax2b.set_ylabel(r'$\Omega$', size = '18')
plt.xlabel(r'Oscillator $i$', size = '18')
# ax2b.set_ylim(-2.5,0)
# ax2b.set_xlim(0,25)
# ax2b.axis("off")

"""
poncaire map
w = 2.8
"""
ax3 = fig.add_subplot(4,2,(2,4))


nspikes = 25   
for s in range(num):
    x_un = spikes_un_w3[s]
    x_sin = spikes_sin_w3[s]
    ax3.plot(x_un[:nspikes], np.zeros(nspikes)+0.25*s+0.15, '|', markersize=18, markeredgewidth=4, markeredgecolor=colors[s])
for s in range(3):
    ax3.plot(x_sin[:nspikes], np.zeros(nspikes)-0.25*s-0.5, '|', markersize=18, markeredgewidth=4, markeredgecolor=mygrey)
for p in range(nspikes):
    plt.axvline(x_sin[p],linestyle = 'dotted', color = 'grey',linewidth=2)

# ax3.set_ylim(-2.5,0)
ax3.set_xlim(0,25)
# ax3.axis("off")

"""
mpv
w = 2.8
"""
ax3b = fig.add_subplot(4,2,6)

plt.plot(np.arange(len(ring_mpv_w3)), ring_mpv_w3, marker = 'o', markersize=3, color = 'k')
plt.xlabel(r'Oscillator $i$', size = '18')
# ax2b.set_ylim(-2.5,0)
# ax2b.set_xlim(0,25)
# ax2b.axis("off")


plt.tight_layout()

today = date.today()
figure_nom = 'FIG4_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()