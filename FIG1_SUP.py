#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:10:22 2022
@author: masoliverm

"""

import matplotlib.pyplot as plt
import numpy as np 
import os
import networkx as nx
from datetime import date

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files/two_pop_chimera'
os.chdir(path)


"""
****************************  DATA   ****************************************
"""

"""
NETWORK drawing: 2 pop
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
twopop chimera
"""
# path_foothills = '/home/maria/cluster_CSG/Documents/Programs/4_FIGURES_ARTICLE_2022'
# path_casa = '/Users/maria/cluster_CSG/Documents/Programs/4_FIGURES_ARTICLE_2022'
# path_csg = '/home/masoliverm/Documents/Programs/4_FIGURES_ARTICLE_2022/'


# os.chdir(path_csg)

K = 0.1
beta = 0.025
ic = 7
nodes = 25
llavor =3
m = 2*nodes
w = 2.6

#path = '/home/maria/cluster_CSG/Documents/project_phase_precession'
#os.chdir(path)


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


"""1. Poincare map"""

unsin = theta
sin = phi

spikes_un = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(unsin)-1):
        if (unsin[i,j] < 0 and unsin[i+1,j] > 0):
            spikes = np.append(spikes,time[i])
    spikes_un.append(spikes)

spikes_sin = []
for j in range(nodes):
    spikes = np.array([])
    for i in range(len(sin)-1):
        if (sin[i,j] < 0 and sin[i+1,j] > 0):
            spikes = np.append(spikes,time[i])
    spikes_sin.append(spikes)


#%%
#from mpl_toolkits.axisartist.axislines import Axes #not working for CSG spyder
import matplotlib.cm as cm
"""
****************************  PLOTS  ****************************************
"""

num = nodes
colors = []
cmap = cm.get_cmap('YlOrRd')
step = (1-0)/25
for i in range(25):
    ind = 0.3+i*step
    colors.append(cmap(ind))

lighteal = '#74a9cf'    
teal = '#3690c0'
darkteal = '#034e7b'
mygrey = '#969696'

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
fig = plt.figure(figsize = (15,5))



"""
131 NETWORK Plotting: 2 pop
"""

ax = fig.add_subplot(131)
options = {"node_size": 400}
optionsS = {
    "node_size": 300,
    "node_color": "white",
    "edgecolors": mygrey,
    "linewidths": 3,
}
optionsD = {
    "node_size": 400,
    "node_color": "white",
    "edgecolors": darkteal,
    "linewidths": 3,
}

nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_1, node_shape = 'd',**optionsD)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=nodelist_2,node_shape = 's',**optionsS)
nx.draw_networkx_nodes(twopop, twopoppos, nodelist=[nodej], node_color=darkteal,node_shape = 'd',**options)

edgesjtp = []
for i in nodelist_1:
    edgesjtp.append((nodej, i))
    
edgesjtp_extra = []
for i in nodelist_2:
    edgesjtp_extra.append((nodej, i))

for i in edgesjtp:
    nx.draw_networkx_edges(twopop, twopoppos, edgelist=[i], width=5 , edge_color=lighteal)
    
for i in edgesjtp_extra:
    nx.draw_networkx_edges(twopop, twopoppos, edgelist=[i], width=1.5 , edge_color=lighteal)

ax.axis("off")




"""
132 Time-series: 2 pop
"""
ax5 = fig.add_subplot(132)
#ax5.set_title('Two Populations: Time-series \n ')

#sax5.text(4490,52,'D',color='k',weight='bold',fontsize = 22)
i = 0
for j in range (9,5,-1):  
    print(j)
    ax5.plot(time, sin[:,j]+10*(i), color=mygrey, linewidth=4)
    ax5.plot(time, unsin[:,j]+10*(i+5), color=mygrey, linewidth=4)
    i = i +1
ax5.plot(time, unsin[:,6]+10*(3+5), color=darkteal, linewidth=4)
ax5.set_xlabel('time $t$')   
ax5.set_xlim(2000,2040)
# ax5.axis("off")

ax3 = fig.add_subplot(133)
nspikes = 800
x_sin = spikes_sin[0]
for p in range(nspikes):
    plt.axvline(x_sin[p],linestyle = 'dotted', color = 'grey',linewidth=2)
i = 0
for s in range(9,5,-1):
    x_un = spikes_un[s]
    ax3.plot(x_un[:nspikes], np.zeros(nspikes)+(i*0.5+3), '|', markersize=18, markeredgewidth=4, markeredgecolor=mygrey)
    x_sin = spikes_sin[s]
    ax3.plot(x_sin[:nspikes], np.zeros(nspikes)+(i*0.5)+0.5, '|', markersize=18, markeredgewidth=4, markeredgecolor=mygrey)
    i = i +1
ax3.plot(x_un[:nspikes], np.zeros(nspikes)+(3*0.5+3), '|', markersize=18, markeredgewidth=4, markeredgecolor=darkteal)    
ax3.plot(time[:len(sin[:,0])], (np.cos(sin[:,0]))/5, markersize=10,  color='#ebbd30' ,linewidth=4)


# ax3.set_ylim(-4,0)
ax3.set_xlim(2000,2040)

# ax3.axis("off")

plt.tight_layout()



today = date.today()


path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)

figure_nom = 'FIG1_supplementary_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()


