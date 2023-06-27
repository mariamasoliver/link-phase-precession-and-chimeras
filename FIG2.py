#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:50:05 2023

@author: maria
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np 
import os
from datetime import date

"""
NETWORK drawing: ring
"""

#nonlinear coupling function
def non_local_coupling(J, j, K):
    coup = np.zeros(J)
    coup[0] = (1+K*np.cos(2*np.pi*abs(j)/J))
    for z in range (1,J):
        coup[z] = coup[z-1] + (1+K*np.cos(2*np.pi*abs(j-z)/J))
    return coup/J

#selected node
nodej = 5
J = 10
K = 0.95

edgesj=[]
for i in range(1,int(J/2)+1):
    edgesj.append((nodej, nodej-i))

for i in range(int(J/2)-1,0,-1):
    edgesj.append((nodej, nodej+i))


coup = non_local_coupling(J, nodej, K)


"""
NETWORK drawing: 2 pop
"""
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
# ti = 1000
# tf = 3000

dtfile = 0.1
# nti = int(ti/dtfile)
# ntf = int(tf/dtfile)

path_files = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files'
os.chdir(path_files)
w = 2.8
data=np.loadtxt('theta_A_0.95_N_500_beta_0.2_ic_chimera_om_2.8_tfinal_10000.txt')
dtfile = 0.1
#we use the same delta we used to compute the mpv


theta_ring = data[:,1:]
time_ring = np.arange(0,len(theta_ring),dtfile)


"""
Loading data 
twopop chimera
"""
path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files/two_pop_chimera'
os.chdir(path)

K = 0.1
beta = 0.025
ic = 7
nodes = 25
llavor =3
m = 2*nodes
w = 2.6

"""Time-series: Removing transients and some extra"""
ti = 1000
tf = 9000
dtfile = 0.1
nti = int(ti/dtfile)
ntf = int(tf/dtfile)
data=np.loadtxt('theta_A_0.1_N_25_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)
data2=np.loadtxt('phi_A_0.1_N_25_beta_0.025_ic_0_omega_'+str(w)+'.txt',skiprows = nti, max_rows = ntf)


theta = data[:,1:]
phi = data2[:,1:]
time = data2[:,0]



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
fig = plt.figure(figsize = (15,10))


"""
231 NETWORK Plotting: 1 pop
"""
myred = '#bb566b'
ax = fig.add_subplot(231)


options1 = {
    "node_size": 300,
    "node_color": "grey",
    "edgecolors": "grey",
    "linewidths": 1,
}
optionsj = {
    "node_size": 300,
    'node_color':myred,
}

nx.draw_networkx_edges(twopop, twopoppos, edgelist=intra_edges, width=0.5, edge_color='grey')
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_1, **options1)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_2, **options1)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=[nodej],**optionsj)
# edges
nx.draw_networkx_edges(
    twopop,
    twopoppos,
    edgelist=inter_edges1,
    width=0.5,
    edge_color='grey',
)
nx.draw_networkx_edges(
    twopop,
    twopoppos,
    edgelist=inter_edges2,
    width=0.5,
    edge_color='grey',
)

for i in range(len(edgesj)):
    nx.draw_networkx_edges(twopop, twopoppos, edgelist=[edgesj[i]], width=coup[i]*5 , edge_color=myred)

plt.axis("off")

"""
234 NETWORK Plotting: 2 pop
"""
lighteal = '#74a9cf'  
darkteal = '#034e7b'
 
ax1 = fig.add_subplot(234)
options = {"node_size": 300}
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_1, node_color='grey',node_shape = 'd',**options)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_2, node_color='grey',node_shape = '^',**options)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=[nodej], node_color=darkteal,node_shape = 'd',**options)


# edges
nx.draw_networkx_edges(
    twopop,
    twopoppos,
    edgelist=inter_edges1,
    width=0.5,
    edge_color='grey',
)
nx.draw_networkx_edges(
    twopop,
    twopoppos,
    edgelist=inter_edges2,
    width=0.5,
    edge_color='grey',
)

 

edgesjtp = []
for i in nodelist_1:
    edgesjtp.append((nodej, i))
    
edgesjtp_extra = []
for i in nodelist_2:
    edgesjtp_extra.append((nodej, i))

nx.draw_networkx_edges(twopop, twopoppos, edgelist=intra_edges, width=0.5, edge_color='grey')
for i in edgesjtp:
    nx.draw_networkx_edges(twopop, twopoppos, edgelist=[i], width=5 , edge_color=darkteal)
    
for i in edgesjtp_extra:
    nx.draw_networkx_edges(twopop, twopoppos, edgelist=[i], width=1.5 , edge_color=darkteal)


ax1.axis("off")

"""
232 SNAPSHOT: 1 pop
"""
nt = 2*10000+500

ax3 = fig.add_subplot(232)
ax3.scatter(np.arange(1,501), theta_ring[nt,:], s =9**2 , color = myred)     
ax3.set_ylabel(r'$\phi$')
# ax2.set_ylim(-3.5,3.5)
ax3.set_xlabel('Nodes $i$') 
     
"""
235 SNAPSHOT: 1 pop
"""
# ti = 5044
ti = 5944

dtplot = 0.1
nti = int(round(ti/dtplot,0))
ntf = int(round(tf/dtplot,0))
time2 = np.arange(ti,tf,dtplot) 
y = np.zeros(6)
for j in range (3):
    y[j] = theta[nti,j]
    y[j+3] = phi[nti,j]
    
ax4 = fig.add_subplot(235)
ax4.set_ylim(-3.5,3.5)
ax4.set_xlabel('Nodes $i$')     
ax4.scatter(np.arange(1,4), y[:3] , marker = 'd', s =12**2 , color = darkteal)     
ax4.scatter(np.arange(4,7), y[3:], marker = '^', s = 12**2, color = darkteal) 
ax4.set_ylabel(r'$\phi$')
ax4.set_xticks([1,2,3,4,5,6])
ax4.set_xticklabels([1,2,3,1,2,3])    

"""
235 
"""
ax5 = fig.add_subplot(233)
ax5.set_ylabel(r'Probability')

"""
236
"""
ax6 = fig.add_subplot(236)
ax6.set_ylabel(r'Probability')
path_fig = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path_fig)
today = date.today()
figure_nom = 'FIG2_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()