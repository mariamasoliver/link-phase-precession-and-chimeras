#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 12:30:07 2023

@author: maria
"""

import matplotlib.pyplot as plt
import numpy as np 
import os
import matplotlib.cm as cm
from datetime import date
import networkx as nx


path = '/Users/maria/Documents/TempCSG/project_spiking_RNN/N_3'
os.chdir(path)


"""
NETWORK construction (as a drawing)
"""

p = 0.05
g = 1.5

np.random.seed(3)
N = 20
A = np.random.normal(0., 1.0, (N,N)) #NxN Matrix with a normal distribution
B = np.random.rand(N,N)<p #NxN Matrix with a uniform distribution. 
                          #It will return true and false, depending on p. 
omega = g*np.multiply(A,B)/np.sqrt(N*p) #Initial weight matrix
positiu = np.where(omega>0)
negatiu = np.where(omega<0)

edges_inh = []
for i in range(1,int(N/2)):
    edges_inh.append((0,i))
    
edges_exc = []
for i in range(int(N/2)+1,N):
    edges_exc.append((10,i))

for i in range(int(N/2)+1,N):
    edges_inh.append((0,i))

for i in range(0,int(N/2)):
    edges_exc.append((10,i))

G_initial = nx.Graph()
G_initial.add_edges_from(edges_exc)
G_l = nx.Graph()
G_l.add_edges_from(edges_inh)
G_initial.add_nodes_from(G_l)

path = '/Users/maria/Documents/TempCSG/project_spiking_RNN/N_3'
os.chdir(path)
vinh = np.loadtxt('vinh_chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt')
vexc = np.loadtxt('vexc_chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt')
fr = np.loadtxt('r_chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt')
output = np.loadtxt('chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt')

nodes = 3

sintheta = output[:,0:nodes] 
costheta = output[:,nodes:2*nodes] 
sinphi = output[:,2*nodes:3*nodes] 
cosphi = output[:,3*nodes:4*nodes] 

theta_o = np.arctan2(sintheta,costheta)
phi_o = np.arctan2(sinphi,cosphi)



step_plot = 2
dtplot = 0.04*step_plot
time_v = np.arange(0,len(vinh)*dtplot,dtplot)
time_r =  np.arange(0,len(fr)*dtplot,dtplot)

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


from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

# VERSION 1: OUTPUT I SUPERVISOR ARE COSPHI COSTHETA 
path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)
fig = plt.figure(figsize = (15,10))

green = '#004D40'
yellow = '#C59606'
darkteal = '#034e7b'

color_edges_e = []
for i in range(len(edges_exc)):
    color_edges_e.append(green)

color_edges_i = []
for i in range(len(edges_inh)):
        color_edges_i.append(yellow)

options_e = {
    "node_size": 350,
    "node_color": "white",
    "edgecolors": green,
    "linewidths": 2,
}

options_i = {
    "node_size": 350,
    "node_color": "white",
    "edgecolors": yellow,
    "linewidths": 2,
}
pos = nx.circular_layout(G_initial) 
edges_width_i = []
for i in range(len(edges_inh)):
    r = np.random.randint(1,3)
    edges_width_i.append(r)
    
edges_width_e = []
for i in range(len(edges_exc)):
    r = np.random.randint(1,3)
    edges_width_e.append(r)


axb = fig.add_subplot(2,3,1)
nx.draw_networkx_nodes(G_initial, pos, nodelist =[10,11,12,13,14,15,16,17,18,19] , node_shape='^',**options_e)
nx.draw_networkx_nodes(G_initial, pos, nodelist =[0,1,2,3,4,5,6,7,8,9] , node_shape='o',**options_i)

# nx.draw_networkx_nodes(G_initial, pos_half_2, node_shape='o',**options)
for i in range(len(edges_inh)):
    nx.draw_networkx_edges(G_initial, pos, edgelist=[edges_inh[i]], width=edges_width_i[i],alpha=1 ,arrows=True, edge_color=color_edges_i[i])
for i in range(len(edges_exc)):
     nx.draw_networkx_edges(G_initial, pos, edgelist=[edges_exc[i]], width=edges_width_e[i], alpha=1, arrows=True,edge_color=color_edges_e[i])
    
plt.axis("off")

#####################
#####################
####################

ax3 = fig.add_subplot(2,3,2)
for i in range(4):
    ax3.plot(time_v, vexc[:,i]+120*i+500, color=green, linewidth=2)
    ax3.plot(time_v, vinh[:,i]+120*i, color=yellow, linewidth=2)
ax3.set_xlabel('time $t$') 
ax3.set_xlim(1000,1500)
# ax3.axis("off")

#####################
#####################
####################

ax = fig.add_subplot(2,3,3)
for i in range(8):
    ax.plot(time_r, fr[:,i]+0.15*i, color='grey', linewidth=2)
ax.set_xlabel('time $t$') 
ax.set_xlim(1000,1500)
# ax.axis("off")

#####################
#####################
####################c
colors = []
cmap = cm.get_cmap('Blues')
step = (1-0)/10
for i in range(10):
    ind = 0.5+i*step
    colors.append(cmap(ind))
    
    

ax2 = fig.add_subplot(2,3,5)
for j in range(3):
    ax2.plot(time_v, costheta[:,j]+3*(j), color=colors[j] ,linewidth=3) 
    ax2.plot(time_v, cosphi[:,j]+3*(j+3), color=colors[j+3] , linewidth=3) 

plt.xlabel(r' time $t$')
ax2.set_xlim(1000,1500)
# ax2.axis("off")

plt.tight_layout()

ax2 = fig.add_subplot(2,3,6)
for j in range(3):
    ax2.plot(time_v, theta_o[:,j]+8*(j), color=colors[j] ,linewidth=3) 
    ax2.plot(time_v, phi_o[:,j]+8*(j+3), color=colors[j+3] , linewidth=3) 

plt.xlabel(r' time $t$')
ax2.set_xlim(1000,1500)
# ax2.axis("off")

plt.tight_layout()

today = date.today()
figure_nom = 'FIG7_SUP_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()

