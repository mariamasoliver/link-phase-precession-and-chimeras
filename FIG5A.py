#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 10:38:55 2023

@author: maria
"""


import matplotlib.pyplot as plt
import numpy as np 
import os
import pandas as pd
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
N = 10

# nodelist_e =[10,11,12,13,14,15,16,17,18,19] #excitatory
# nodelist_i =[0,1,2,3,4,5,6,7,8,9] #inhibitory

nodelist_e =[5,6,7,8,9] #excitatory
nodelist_i =[0,1,2,3,4] #inhibitory


edges_inh = []
for j in range(0,N):
    for i in nodelist_i:
        if(j != i):
            edges_inh.append((j,i))
            
edges_exc = []
for j in range(0,N):
    for i in nodelist_e:
        if(j != i):
            edges_exc.append((j,i))


G_initial = nx.Graph()
G_initial.add_edges_from(edges_exc)
G_initial.add_edges_from(edges_inh)
# G_l = nx.Graph()
# G_l.add_edges_from(edges_inh)
# G_initial.add_nodes_from(G_l)

path = '/Users/maria/Documents/TempCSG/project_spiking_RNN/N_3'
os.chdir(path)
vinh = np.loadtxt('vinh_chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt')
vexc = np.loadtxt('vexc_chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt')


step_plot = 2
dtplot = 0.04*step_plot
time_v = np.arange(0,len(vinh)*dtplot,dtplot)


#############
#############


#Read the firing times & the sincronized output
n_tspike_step_1 = np.load('n_tspike_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.npy', allow_pickle=True)
output = np.loadtxt('chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt') 

#define some parameters
tmax = 10000
step = 2
dt = 0.04
nmax = tmax/dt
time_vec = np.arange(0,tmax+0.04,dt*step)
N = 10000
sintheta = output[:,0] 

# 1. Compute the hilbert transform of the synchronized activity to extract the phase, and extract the phases
# that correspond to the firing times
from scipy.signal import hilbert
rphase = np.real(hilbert(sintheta))
iphase = np.imag(hilbert(sintheta)) #form the analytical signal
ht_phase = np.arctan2(iphase,rphase)#inst phase

# convert spiking times to the other sampling rate
#1. Convert spiking times times
n_tspike=[]

for i in range(N):
     n_tspike.append((n_tspike_step_1[i]/step).astype(int))
     
# 2. obtaining the corresponding phase of the firing times
phase_firing_ntimes = []
firing_ntimes = []
av_phase = np.zeros(N)

for i in range (N):
    phase_firing_ntimes.append(ht_phase[n_tspike[i]])
    firing_ntimes.append(n_tspike[i]) #just reformating it as list of arrays instead of an array of arrays-easier for dataframes

p_min = np.zeros(N)
p_max = np.zeros(N)

for n in range(0,N):
    if (len(phase_firing_ntimes[n])==0):
        p_min[n] = 0
        p_max[n] = 0
    else:
        p_min[n] = phase_firing_ntimes[n].min()
        p_max[n] = phase_firing_ntimes[n].max()
        
## 1. Sort neurons in function of their phase preference
### We can clearly see two diferent phase distributions: modal (or bidmodal) and uniformal
### We are going to define the neuron classifier index. It is basically a re-definition of the order parameter. 
# It quantifies how much the phase distribution for a given neuron is spread (small neuron classifier index) 
# or not (large neuron classifier index). 

def NCI(phase):
    """Returns the neuron classifier index for a 1xT numpy array.
    Parameters:
    -----------
    phase: 1xT numpy array (T: number of spikes).
               Each row is the time the neuron fired
    Returns:
    --------
    r: Neuron classifier index (Scalar number at each time)
    """
    scos = np.cos(phase).sum()
    ssin = np.sin(phase).sum()
    r = np.sqrt((scos*scos + ssin*ssin)) / (len(phase))
    return (np.round(r,6))    

### 1.1 Excitatory neurons and inhibitory

Nsub = 5000
NCIvecE = np.zeros(Nsub)
for n in range(Nsub):
    phase = phase_firing_ntimes[n]
    if (len(phase) == 0):
        NCIvecE[n] = np.nan
        print(n)
    else:
        NCIvecE[n] = NCI(phase)



NCIpdE = pd.DataFrame(NCIvecE)
print(NCIpdE.shape)
NCIpdENoNans = NCIpdE.dropna()
print(NCIpdENoNans.shape)

Nsub = 5000
NCIvecI = np.zeros(Nsub)
for n in range(Nsub):
    phase = phase_firing_ntimes[n+Nsub]
    if (len(phase) == 0):
        NCIvecI[n] = np.nan
        print(n)
    else:
        NCIvecI[n] = NCI(phase)



NCIdic = {}
NCIdic['NCI Exc'] = NCIvecE
NCIdic['NCI Inh'] = NCIvecI
NCIpd = pd.DataFrame(NCIdic)

NCIpdNoNans = NCIpd.dropna()

#Sorting neurons by using the NCI
neurons_data = pd.DataFrame(NCIpd)

neurons_data['Exc Class'] = neurons_data['NCI Exc']<=0.6
neurons_data['Inh Class'] = neurons_data['NCI Inh']<=0.6

Esorted = neurons_data.sort_values(["NCI Exc"], ascending = True)
Isorted = neurons_data.sort_values(["NCI Inh"], ascending = True)

Esorted = Esorted.drop(['NCI Inh', 'Inh Class'], axis=1)
Isorted = Isorted.drop(['NCI Exc', 'Exc Class'], axis=1)

print(Esorted.head())
print(Isorted.head())

#save the phases on a dataframe, each row a column
phases = pd.DataFrame(phase_firing_ntimes).T
firing_ntimes_df = pd.DataFrame(firing_ntimes).T

#save the indices depending on their class type and phase precession 
Eind_pp = np.array(Esorted[Esorted['Exc Class']==True]['Exc Class'].index)
Eind_nopp = np.array(Esorted[Esorted['Exc Class']==False]['Exc Class'].index)

Iind_pp = np.array(Isorted[Isorted['Inh Class']==True]['Inh Class'].index)
Iind_nopp = np.array(Isorted[Isorted['Inh Class']==False]['Inh Class'].index)




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
fig = plt.figure(figsize = (25,5))

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
    "node_color": green,
    "edgecolors": green,
    "linewidths": 2,
}

options_i = {
    "node_size": 350,
    "node_color": yellow,
    "edgecolors": yellow,
    "linewidths": 2,
}
pos = nx.circular_layout(G_initial) 
edges_width_i = []
seed = 28
np.random.seed(seed = seed)
for i in range(len(edges_inh)):
    r = np.random.randint(1,3)
    edges_width_i.append(r+1)
    
edges_width_e = []
for i in range(len(edges_exc)):
    r = np.random.randint(1,3)
    edges_width_e.append(r)


axb = fig.add_subplot(1,5,1)
nx.draw_networkx_nodes(G_initial, pos, nodelist =[5,6,7,8,9] , node_shape='^',**options_e)
nx.draw_networkx_nodes(G_initial, pos, nodelist =[0,1,2,3,4] , node_shape='o',**options_i)

# nx.draw_networkx_nodes(G_initial, pos_half_2, node_shape='o',**options)
for i in range(len(edges_inh)):
    nx.draw_networkx_edges(G_initial, pos, edgelist=[edges_inh[i]], width=edges_width_i[i],alpha=1 ,arrows=False, edge_color=color_edges_i[i])
for i in range(len(edges_exc)):
     nx.draw_networkx_edges(G_initial, pos, edgelist=[edges_exc[i]], width=edges_width_e[i], alpha=1, arrows=False,edge_color=color_edges_e[i])
    
plt.axis("off")

#####################
#####################
####################

ax3 = fig.add_subplot(1,5,2)
for i in range(4):
    ax3.plot(time_v, vexc[:,i]+120*i+500, color=green, linewidth=2)
    ax3.plot(time_v, vinh[:,i]+120*i, color=yellow, linewidth=2)
ax3.set_xlabel('time $t$') 
ax3.set_xlim(1000,1500)
# ax3.axis("off")

#####################
#####################
####################c
ax4 = fig.add_subplot(1,5,3)
ax4 = fig.add_subplot(1,5,4)   
ax4 = fig.add_subplot(1,5,5)   


plt.tight_layout()

today = date.today()
figure_nom = 'FIG5A_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()

