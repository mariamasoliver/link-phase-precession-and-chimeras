#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:37:22 2023
@author: masoliverm

"""

import matplotlib.pyplot as plt
import numpy as np 
import os
import networkx as nx
from datetime import date
from scipy.signal import find_peaks

"""
****************************  DATA   ****************************************
"""


"""
NETWORK drawing: 1 pop
"""

#nonlinear coupling function
def non_local_coupling(J, j, K):
    coup = np.zeros(J)
    for z in range (J):
        coup[z] = (1+K*np.cos(2*np.pi*abs(j-z)/J))
    return coup

#selected node
nodej = 5
J = 10
K = 0.95

edgesj=[]
for i in range(1,int(J/2)+1):
    edgesj.append((nodej, nodej-i))

for i in range(int(J/2)-1,0,-1):
    edgesj.append((nodej, nodej+i))


coup = non_local_coupling(J, 0, K)
#%%
twopop = nx.complete_graph(J)
twopoppos = nx.circular_layout(twopop)  # positions for all nodes

# nodes
nodelist_1 = [3,4,5,6,7]
nodelist_2 = [2,1,0,9,8]

intra_edges=[]
for i in range(len(nodelist_1)):
    for j in range(len(nodelist_2)):
        intra_edges.append((nodelist_1[i], nodelist_2[j]))
        
        
inter_edges1=[]
for i in nodelist_1:
    for j in nodelist_1:
        inter_edges1.append((i, j))

inter_edges2=[]
for i in nodelist_2:
    for j in nodelist_2:
        inter_edges2.append((i, j))

     

"""
Loading data 
onepop chimera
"""
ti = 1000
tf = 3000
dtfile = 0.1
nti = int(ti/dtfile)
ntf = int(tf/dtfile)

path_files = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files'
os.chdir(path_files)

data=np.loadtxt('theta_A_0.95_N_500_beta_0.2_ic_chimera_om_1.0_tfinal_10000.txt',skiprows = nti, max_rows = ntf)
dtfile = 0.1
#we use the same delta we used to compute the mpv
deltat = 1000
ndeltat = int(deltat/dtfile)
theta_ring = data[ndeltat:2*ndeltat,1:]
time_ring = data[ndeltat:2*ndeltat,0]

"""Poincare map"""

#need to know which are the sincro/unsincro neurons
#we'll use the mpv. The first 10 nodes with the lowest mpv: sincro pop
#the last 10 nodes with the highest mpv: unsincro pop

ring_mpv = np.loadtxt('mpv_omega_1.0_N_500.txt') 

nodes = 10
mpv_max = np.max(ring_mpv)
mpv_min = np.min(ring_mpv)

mpv_max_node = np.where(ring_mpv==mpv_max)
mpv_min_node = np.where(ring_mpv==mpv_min)

unsin = theta_ring[:,mpv_max_node[0][0]:mpv_max_node[0][0]+nodes]
sin = theta_ring[:,mpv_min_node[0]]


spikes_un = []
for j in range(nodes):
    spikes = np.array([])
    peaks_t = find_peaks(unsin[:,j],height=2.5)
    spikes_un.append(time_ring[peaks_t[0]])

spikes_sin = []
for j in range(nodes):
    spikes = np.array([])
    peaks_t = find_peaks(sin[:,j],height=2.5)
    spikes_sin.append(time_ring[peaks_t[0]])


#%%
#from mpl_toolkits.axisartist.axislines import Axes #not working for CSG spyder
import matplotlib.cm as cm
"""
****************************  PLOTS  - OPTION 2 ****************************************
"""
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
fig = plt.figure(figsize = (20,5))

"""
152 NETWORK Plotting: 1 pop
"""
myred = '#bb566b'
mylightred = '#D699A6'
mygrey = '#969696'
ax = fig.add_subplot(141)


options1 = {
    "node_size": 400,
    "node_color": "white",
    "edgecolors": mygrey,
    "linewidths": 3,
}
optionsj = {
    "node_size": 400,
    'node_color':myred,
    "edgecolors": myred,
    "linewidths": 3,
}

nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_1, **options1)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_2, **options1)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=[nodej],**optionsj)
# edges


for i in range(len(edgesj)):
    nx.draw_networkx_edges(twopop, twopoppos, edgelist=[edgesj[i]], width=coup[i]*3 , edge_color=mylightred)

plt.axis("off")





"""
15(3,4) Time-series: 1 pop
"""

ax3 = fig.add_subplot(1,4,2)
#ax5.set_title('Two Populations: Time-series \n ')

#sax5.text(4490,52,'D',color='k',weight='bold',fontsize = 22)
for j in range (0,4):  
    ax3.plot(time_ring, unsin[:,j]+8*(j+4), color=mygrey, linewidth=4)
    ax3.plot(time_ring, sin[:,j]+8*(j), color=mygrey, linewidth=4)
    
ax3.plot(time_ring, unsin[:,3]+8*(3+4), color=myred, linewidth=4)
ax3.plot(time_ring, sin[:,0], color=mygrey, linewidth=4)

ax3.set_xlabel('time $t$') 
# ax3.set_yticklabels([])
# ax3.set_yticks([])    
ax3.set_xlim(2100,2125)

# ax3.axis("off")


"""
155 Poincaire-like map 1 pop
"""
ax5 = fig.add_subplot(1,4,4)

num = nodes
r = np.random.uniform(size=num)

x_sin = spikes_sin[0]
nspikes = 25 
for p in range(nspikes):
    plt.axvline(x_sin[p],linestyle = 'dotted', color = 'grey',linewidth=2)

 
for s in range(0,10):
    x_un = spikes_un[s]
    x_sin = spikes_sin[s]
    ax5.plot(x_un[:nspikes], np.zeros(nspikes)+0.35*s, '|', markersize=22, markeredgewidth=5, markeredgecolor=mygrey)
    # ax5.plot(x_sin[:nspikes], np.zeros(nspikes)-0.35*s-0.35, '|', markersize=22, markeredgewidth=5, markeredgecolor=mygrey)
x_un = spikes_un[s]
ax5.plot(x_un[:nspikes], np.zeros(nspikes)+0.35*s, '|', markersize=22, markeredgewidth=5, markeredgecolor=myred)



# ax5.plot(time_ring, (np.cos(sin[:,0]+np.pi))/8 -1.8, markersize=10,  color='#ebbd30' ,linewidth=4)
# ax5.set_ylim(-3,0)
ax5.set_xlim(2100,2200)
# ax5.axis("off")



plt.tight_layout()
today = date.today()

path_fig = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path_fig)

figure_nom = 'FIG3_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()

