#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:10:19 2023

@author: maria
"""
import numpy as np
import os
import pandas as pd
path = '/Users/maria/Documents/TempCSG/project_spiking_RNN/N_3'
os.chdir(path)

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
#Let's check the time series
#showing phase precession
ntoplot = 5
import seaborn as sns
import matplotlib.pyplot as plt
colors = sns.color_palette("dark:salmon_r", ntoplot)
Ninh = 5000

i = 1
fig = plt.figure(figsize = (20,2*ntoplot))
for n in Iind_pp[0:ntoplot]:
    ax = fig.add_subplot(ntoplot,1,i)
    ax.plot(time_vec[n_tspike[n+Ninh]], phase_firing_ntimes[n+Ninh], '.', markersize=15, color=colors[0])
    # ax.plot(time_vec[n_tspike[n]], phase_firing_ntimes[n], '.', markersize=15, color=colors[0])
    ax.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =2)
    i = i+1

    ax.set_xlim(400,3000)
    
#%%
#Let's check the time series
#not showing phase precession
ntoplot = 5
import seaborn as sns
import matplotlib.pyplot as plt
colors = sns.color_palette("dark:salmon_r", ntoplot)
Ninh = 5000

i = 1
fig = plt.figure(figsize = (20,2*ntoplot))
for n in Iind_nopp[-ntoplot:]:
    ax = fig.add_subplot(ntoplot,1,i)
    ax.plot(time_vec[n_tspike[n+Ninh]], phase_firing_ntimes[n+Ninh], '.', markersize=15, color=colors[0])
    # ax.plot(time_vec[n_tspike[n]], phase_firing_ntimes[n], '.', markersize=15, color=colors[0])
    ax.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =2)
    i = i+1

    ax.set_xlim(400,3000)
    
#%%

#We are going to chose 10 neurons to show phase precession
#let's define a trial using the plot of the theta oscillation
#find where the theta has its maxs. Every max will be a trial
#plot it between 0 and 3pi as wilten suggested

from scipy.signal import find_peaks
peaks_syn = []
peaks_t = find_peaks(ht_phase,height=1)
itrials = peaks_t[0]
ntrials = len(itrials)
ttrials = time_vec[peaks_t[0]]


#sanity check
fig = plt.figure(figsize = (20,2.5))
ax = fig.add_subplot(1,1,1)
ax.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =2)
ax.plot(ttrials, ht_phase[peaks_t[0]], 'x', mew = 4, color ='grey')

ax.set_ylabel(r'$phase$')
ax.set_xlabel('Time (ms)')
ax.set_xlim(1000,3000)


#%%
colors_per_trial = ['#fdae6b', '#fd8d3c', '#f16913','#d94801','#8c2d04']
#split firing times into trials
#for 10 neurons showing phase precession
for n in Iind_pp[:ntoplot]:
    ntspike_by_trial = []
    for i in range(len(itrials)-1):
        ntimes = np.array([])
        for j in n_tspike[n+Ninh]:
            if ((j > itrials[i]) & (j < itrials[i+1])):
                ntimes = np.append(ntimes,j)
        ntspike_by_trial.append(ntimes.astype(int))

    ntspike_by_trial_shifted = []
    for i in range(len(itrials)-1):
        ntspike_by_trial_shifted.append(ntspike_by_trial[i][:]-itrials[i])


    phase_by_trial = []  
    #get the phase for those times
    for j in range(0,len(itrials)-1):
        p = np.array([])
        for i in ntspike_by_trial[j]:
            p = np.append(p,ht_phase[i])
        phase_by_trial.append(p)
#PLOT JUST FOR THE FIRST X TRIALS
    ntrials_to_plot = 5
    fig = plt.figure(figsize = (10,5))
    ax = fig.add_subplot(111)
    for j in range(ntrials_to_plot):
        ax.plot(ntspike_by_trial_shifted[j][0],phase_by_trial[j][0],'o',color=colors_per_trial[j])
        #ax.plot(ntspike_by_trial_shifted[j][0],phase_by_trial[j][0]+3*np.pi,'o',color=colors_per_trial[j])
        #ax.plot(phase_by_trial[j]+np.pi,'o',color=yellow)
    ax.set_ylim(-np.pi,np.pi)
    ax.set_xlabel('Shifted time (s)')
    ax.set_ylabel('Phase')


#%%

#We are going to chose 10 neurons NOT SHOWING phase precession
#let's define a trial using the plot of the theta oscillation
#find where the theta has its maxs. Every max will be a trial
#plot it between 0 and 3pi as wilten suggested


#split firing times into trials
#for 10 neurons NOT showing phase precession
for n in Iind_nopp[-ntoplot:]:
    ntspike_by_trial = []
    for i in range(len(itrials)-1):
        ntimes = np.array([])
        for j in n_tspike[n+Ninh]:
            if ((j > itrials[i]) & (j < itrials[i+1])):
                ntimes = np.append(ntimes,j)
        ntspike_by_trial.append(ntimes.astype(int))

    ntspike_by_trial_shifted = []
    for i in range(len(itrials)-1):
        ntspike_by_trial_shifted.append(ntspike_by_trial[i][:]-itrials[i])


    phase_by_trial = []  
    #get the phase for those times
    for j in range(0,len(itrials)-1):
        p = np.array([])
        for i in ntspike_by_trial[j]:
            p = np.append(p,ht_phase[i])
        phase_by_trial.append(p)
    
    ntrials_to_plot = 5
    fig = plt.figure(figsize = (10,5))
    ax = fig.add_subplot(111)
    for j in range(ntrials_to_plot):
        ax.plot(ntspike_by_trial_shifted[j][0],phase_by_trial[j][0],'o',color = colors_per_trial[j])
        #ax.plot(tspike_by_trial_shifted[j],phase_by_trial[j]+3*np.pi,'o',color=yellow)
        #ax.plot(phase_by_trial[j]+np.pi,'o',color=yellow)
    ax.set_ylim(-np.pi,np.pi)
    ax.set_xlabel('Shifted time (s)')
    ax.set_ylabel('Phase')




#%%##########################
###########################
###########################

green = '#004D40'
colors = ['#F4B028','#C59606', '#D87E0A']
yellow = ['#fe9929', '#ec7014', '#cc4c02']


params = {
    'axes.labelsize': 18,
        'xtick.labelsize': 18,
            'ytick.labelsize': 18,
                'legend.fontsize': 18,
                'axes.titlesize':18,
                    'text.usetex': False,
                        #'font': 'Helvetica',
                        'mathtext.bf': 'helvetica:bold',

                    }

plt.rcParams.update(params)
fig = plt.figure(figsize = (10,3))
ax = fig.add_subplot(321)
Ninh = 5000
ppi = 2
ax.plot(time_vec[n_tspike[Iind_pp[ppi]+Ninh]], np.zeros(len(time_vec[n_tspike[Iind_pp[ppi]+Ninh]]))+5, '|', markersize=22, mew = 3,color=colors[0])
ax.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =3)
ax.plot(time_vec[n_tspike[Iind_pp[ppi]+Ninh]], phase_firing_ntimes[Iind_pp[ppi]+Ninh], '.', markersize=18, color=colors[0])


ax.set_xlim(980 ,1380)
ax.set_xlabel(' ')
ax.set_xticks([])
plt.axis("off")
ax2 = fig.add_subplot(323)
Ninh = 5000
ppi = 3
ax2.plot(time_vec[n_tspike[Iind_pp[ppi]+Ninh]], np.zeros(len(time_vec[n_tspike[Iind_pp[ppi]+Ninh]]))+5, '|', markersize=22, mew = 3,color=colors[1])
ax2.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =3)
ax2.plot(time_vec[n_tspike[Iind_pp[ppi]+Ninh]], phase_firing_ntimes[Iind_pp[ppi]+Ninh], '.', markersize=18, color=colors[1])

ax2.set_xlim(980,1380)
ax2.set_xlabel(' ')
ax2.set_xticks([])
plt.axis("off")

ax4 = fig.add_subplot(325)
Ninh = 5000
ppi = 1
ax4.plot(time_vec[n_tspike[Iind_pp[ppi]+Ninh]], np.zeros(len(time_vec[n_tspike[Iind_pp[ppi]+Ninh]]))+5, '|', markersize=22, mew = 3,color=colors[2])
ax4.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =3)
ax4.plot(time_vec[n_tspike[Iind_pp[ppi]+Ninh]], phase_firing_ntimes[Iind_pp[ppi]+Ninh], '.', markersize=18, color=colors[2])
ax4.set_xlim(980,1380)
plt.axis("off")

#####################

ax5 = fig.add_subplot(322)
Ninh = 5000
ppi = -2
ax5.plot(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]], np.zeros(len(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]]))+5, '|', markersize=22, mew = 3,color=yellow[0])
ax5.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =3)
ax5.plot(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]], phase_firing_ntimes[Iind_nopp[ppi]+Ninh], '.', markersize=18, color=yellow[0])


ax5.set_xlim(980 ,1380)
ax5.set_xlabel(' ')
ax5.set_xticks([])
plt.axis("off")

ax6 = fig.add_subplot(324)
Ninh = 5000
ppi = -3
ax6.plot(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]], np.zeros(len(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]]))+5, '|', markersize=22, mew = 3,color=yellow[1])
ax6.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =3)
ax6.plot(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]], phase_firing_ntimes[Iind_nopp[ppi]+Ninh], '.', markersize=18, color=yellow[1])

ax6.set_xlim(980,1380)
ax6.set_xlabel(' ')
ax6.set_xticks([])
plt.axis("off")

ax7 = fig.add_subplot(326)
Ninh = 5000
ppi = -1
ax7.plot(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]], np.zeros(len(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]]))+5, '|', markersize=22, mew = 3,color=yellow[2])
ax7.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =3)
ax7.plot(time_vec[n_tspike[Iind_nopp[ppi]+Ninh]], phase_firing_ntimes[Iind_nopp[ppi]+Ninh], '.', markersize=18, color=yellow[2])

ax7.set_xlim(980,1380)

plt.axis("off")

from datetime import date
path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)
today = date.today()
figure_nom = 'FIG5B_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()

#%%
path = '/Users/maria/Documents/TempCSG/project_spiking_RNN/N_3'
os.chdir(path)
vinh = np.loadtxt('vinh_chimera_Q_1400_G_15000_N_10000_T_10000_bias_2000_ic_30_Tt_10000_dl_od.txt')

import matplotlib.pyplot as plt
from scipy.signal import find_peaks
peaks_vinh = find_peaks(vinh[:,0],height=20)


step_plot = 2
dtplot = 0.04*step_plot
time_v = np.arange(0,len(vinh)*dtplot,dtplot)
t_peaks_vinh = time_v[peaks_vinh[0]]


yellow = '#C59606'
params = {
    'axes.labelsize': 18,
        'xtick.labelsize': 18,
            'ytick.labelsize': 18,
                'legend.fontsize': 18,
                'axes.titlesize':18,
                    'text.usetex': False,
                        #'font': 'Helvetica',
                        'mathtext.bf': 'helvetica:bold',

                    }

plt.rcParams.update(params)
fig = plt.figure(figsize = (5,3))
ax = fig.add_subplot(3,1,1)
for i in range(len(peaks_vinh[0][:100])):
    ax.axvline(time_vec[peaks_vinh[0][i]], color = 'grey', ls = '--', lw = '2')
ax.set_xlim(980,1380)
plt.axis("off")
ax.plot(time_v, vinh[:,0], color=yellow, linewidth=3)

ax.set_xlim(980,1380)
plt.axis("off")

ax2 = fig.add_subplot(3,1,2)
for i in range(len(peaks_vinh[0][:100])):
    ax2.axvline(time_vec[peaks_vinh[0][i]], color = 'grey',ls = '--', lw = '2')

plt.axis("off")
ax2.plot(t_peaks_vinh, np.zeros(len(vinh[peaks_vinh[0],0]))+5, '|', markersize=22, mew = 3,color=yellow)
ax2.plot(time_vec[:len(ht_phase)], ht_phase, color='k', linewidth =3)
ax2.plot(time_vec[peaks_vinh[0][:100]], ht_phase[peaks_vinh[0][:100]], '.', markersize=18, color=yellow)



ax2.set_xlabel('time $t$') 
ax2.set_xlim(980,1380)
plt.axis("off")
# ax.plot(t_peaks_vinh, vinh[peaks_vinh[0],0], 'x', mew = 4, color ='grey')

ax3 = fig.add_subplot(3,1,3)
for i in range(len(peaks_vinh[0][:100])):
    ax3.axvline(time_vec[peaks_vinh[0][i]], color = 'grey',ls = '--', lw = '2')
ax3.set_xlim(980,1380)
plt.axis("off")
from datetime import date
path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)
today = date.today()
figure_nom = 'FIG5C_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()

#%%





