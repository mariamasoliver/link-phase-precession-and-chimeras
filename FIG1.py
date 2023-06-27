#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:29:08 2023

@author: maria
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from datetime import date

params = {
    'axes.labelsize': 22,
        'xtick.labelsize': 22,
            'ytick.labelsize': 22,
                'legend.fontsize': 22,
                     'axes.titlesize':20,
                    'text.usetex': False,
                        #'font': 'Helvetica',
                        'mathtext.bf': 'helvetica:bold',

                    }

plt.rcParams.update(params)
dt = 0.1
nt = 660
time = np.arange(-12,nt-12)


colors = ['#66c2a5','#fc8d62','#8da0cb', '#80b1d3']
freq = 1/60
sup = np.sin(2*time*np.pi*freq)
peaks, _ = find_peaks(sup, height=1)
spikes = len(peaks[1:])


"""
Figure Phase precession 1 neuron
"""
fig = plt.figure(figsize = (10,4))
ax = fig.add_subplot(2,1,1)
shift = 1
#ax.text(550,1.8,'Pyramidal neuron',color=colors[0],horizontalalignment='center',verticalalignment='center',fontsize = 18)
#ax.text(550,0,'Theta oscillation',color='#ebbd30',horizontalalignment='center',verticalalignment='center',fontsize = 18)
ax.set_title('Phase precession', fontsize = 14, pad=20)
ax.arrow(-12, -2-shift, 660, 0, head_width=0.3, head_length=8, fc='grey', ec='grey',lw=3)
ax.plot(time, sup-shift, color='#ebbd30', linewidth=3)
ax.text(725,1.8,'Pyramidal neuron',color=colors[0],horizontalalignment='center',verticalalignment='center',fontsize = 12)
ax.text(725,0.3,'Theta oscillation',color='#ebbd30',horizontalalignment='center',verticalalignment='center',fontsize = 12)

ax.text(725,-1.3,'Interneuron',color='#80b1d3',horizontalalignment='center',verticalalignment='center',fontsize = 12)


for p in peaks:
    plt.vlines(x=p-12,ymin=-0.5-shift,ymax= 3.5, linestyle = 'dashed', color = 'k', alpha = 0.4)
precess = np.arange(8,8*len(peaks),8)
ax.plot(peaks[1:]-precess[:]-12, np.zeros(spikes)+1.8, '|', markersize=10, markeredgewidth=3, markeredgecolor=colors[0])
#for s in range(1,4):
    #ax.plot(peaks[1+s:-s]-precess[s:-s]-12+4*s, np.zeros(spikes-2*s)+1.8, '|', markersize=10, markeredgewidth=2, markeredgecolor=colors[0])

ax.plot(peaks-12, np.zeros(len(peaks))-1.3-1, '|', markersize=10, markeredgewidth=3, markeredgecolor='#80b1d3')



ax.set_ylim(-2-shift-1,3)
ax.axis("off")

"""
Figure Phase precession 3 neurons
"""


ax2 = fig.add_subplot(2,1,2)
ax2.set_title('Theta sequences', fontsize = 14, pad=20)

ax2.text(725,1.9,'Pyramidal neuron 1',color=colors[0],horizontalalignment='center',verticalalignment='center',fontsize = 12)
ax2.text(725,1.1,'Pyramidal neuron 2',color=colors[1],horizontalalignment='center',verticalalignment='center',fontsize = 12)
ax2.text(725,0.3,'Pyramidal neuron 3',color=colors[2],horizontalalignment='center',verticalalignment='center',fontsize = 12)

ax2.text(725,-0.7,'Theta oscillation',color='#ebbd30',horizontalalignment='center',verticalalignment='center',fontsize = 12)
ax2.plot(time, sup-shift, color='#ebbd30', linewidth=3)
for p in peaks:
    plt.vlines(x=p-12,ymin=-0.5-shift,ymax= 3.5, linestyle = 'dashed', color = 'k', alpha = 0.4)
precess = np.arange(8,8*len(peaks),8)

ax2.plot(peaks[1:]-precess[:]-12, np.zeros(spikes)+2.8, '|', markersize=10, markeredgewidth=3, markeredgecolor=colors[0])
#for s in range(1,4):
    #ax2.plot(peaks[1+s:-s]-precess[s:-s]-12+4*s, np.zeros(spikes-2*s)+1.8, '|', markersize=10, markeredgewidth=2, markeredgecolor=colors[0])

ax2.plot(peaks[2:]-precess[:-1]-12, np.zeros(len(precess[:-1]))+2.8-0.8, '|', markersize=10, markeredgewidth=3, markeredgecolor=colors[1])
#ax2.plot(peaks[2+1:]-precess[1:-1]-12+4*1, np.zeros(len(precess[1:-1]))+1.8-0.3, '|', markersize=10, markeredgewidth=2, markeredgecolor=colors[1])
#for s in range(2,4):
    #ax2.plot(peaks[2+s:-s+1]-precess[s:-s]-12+4*s, np.zeros(len(precess[s:-s]))+1.8-0.3, '|', markersize=10, markeredgewidth=2, markeredgecolor=colors[1])

ax2.text(725,-1.5,'Interneurons',color='#80b1d3',horizontalalignment='center',verticalalignment='center',fontsize = 12)
for i in range(3):
    ax2.plot(peaks-12-5*i+5, np.zeros(len(peaks))-1.3-1, '|', markersize=10, markeredgewidth=3, markeredgecolor=colors[3])



ax2.plot(peaks[3:]-precess[:-2]-12, np.zeros(len(precess[2:]))+2.8-1.6, '|', markersize=10, markeredgewidth=3, markeredgecolor=colors[2])
#ax2.plot(peaks[3+1:]-precess[1:-2]-12+4*1, np.zeros(len(precess[1+2:]))+1.8-0.9, '|', markersize=10, markeredgewidth=2, markeredgecolor=colors[2])
#ax2.plot(peaks[3+2:]-precess[2:-2]-12+4*2, np.zeros(len(precess[1+3:]))+1.8-0.9, '|', markersize=10, markeredgewidth=2, markeredgecolor=colors[2])
#ax2.plot(peaks[3+3:-1]-precess[3:-3]-12+4*3, np.zeros(len(precess[1+4:-1]))+1.8-0.9, '|', markersize=10, markeredgewidth=2, markeredgecolor=colors[2])
ax2.arrow(-12, -2-shift, 660, 0, head_width=0.3, head_length=8, fc='grey', ec='grey', lw = 3)
ax2.set_ylim(-2-shift-1,3.5)
ax2.axis("off")
plt.tight_layout()
today = date.today()
figure_nom = 'FIG1_pp_drawing'+str(today)+'.svg'
plt.savefig(figure_nom, bbox_inches='tight', format='svg',dpi = 400)
