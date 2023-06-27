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

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files/two_pop_chimera'
os.chdir(path)

"""RATIO MPV"""
onepopratio = np.loadtxt('meanratio_vs_w_N_500_coarse_long.txt')
twopopratio = np.loadtxt('ratio_mpv_A_0.1_delta_1000.0_ic_0_TC_1_N_25_two_pop_coarse.txt')



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
fig = plt.figure(figsize = (12,3.25))

ax1 = fig.add_subplot(1,2,1)

lighteal = '#74a9cf'  
darkteal = '#034e7b'
myred = '#bb566b'

ax1.plot(onepopratio[:,0], onepopratio[:,1], marker = 'o', lw = '3.5', markersize=10, color = myred)
plt.fill_between(onepopratio[:,0], onepopratio[:,1] - onepopratio[:,2], onepopratio[:,1] + onepopratio[:,2],color=myred, alpha=0.4)
ax1.axhline(y = 1,color = 'grey',linestyle = '--', lw = '2')

ax1.axhline(y = 0.88,color = 'grey',linestyle = '--', lw = '2')
plt.xlabel(r'Intrinsic frequency $\rho$')
ax1.set_ylim(0.25,1.1)


ax2 = fig.add_subplot(1,2,2)
ax2.axhline(y = 0.88,color = 'grey',linestyle = '--', lw = '2')
ax2.plot(np.arange(1,11), twopopratio, marker = 'o', lw = '3.5',markersize=10, color = darkteal)
ax1.set_ylabel(r'$\Omega_{ratio}$')

ax2.set_yticklabels([])
plt.xlabel(r'Intrinsic frequency $\rho$')
ax2.axhline(y = 1,color = 'grey',linestyle = '--', lw = '2')
ax2.set_ylim(0.25,1.1)

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022'
os.chdir(path)

plt.tight_layout()
today = date.today()
figure_nom = 'FIG4_sup_'+str(today)+'.svg'
plt.tight_layout()
plt.savefig(figure_nom, bbox_inches='tight',format='svg', dpi=300)
plt.show()
