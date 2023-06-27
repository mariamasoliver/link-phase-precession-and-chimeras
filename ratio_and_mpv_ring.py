# -*- coding: utf-8 -*-
"""
Note: we compute the ratio between the sincronized and the unsyncronized pop
for the ring set-up, which as we know shows a drift. Meaning that for each
delta t the population that is sincronized changes. Here we do it for different
delta t and average.
"""

import numpy as np
import os
from scipy.signal import find_peaks

# params = {
#     'axes.labelsize': 14,
#         'xtick.labelsize': 14,
#             'ytick.labelsize': 14,
#                 'legend.fontsize': 14,
                    
#                     'text.usetex': False,
#                         #'font': 'Helvetica',
#                         'mathtext.bf': 'helvetica:bold',

#                     }


A = 0.95
beta = 0.2
N = 100
omega = 3.2
dt = 0.1

path = '/Users/maria/Documents/TempCSG/4_FIGURES_ARTICLE_2022/files'


w_vec_ring = np.array([1.0, 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.5 ,
        1.55, 1.6 , 1.65, 1.7 , 1.75, 1.8 , 1.85, 1.9 , 1.95, 2.0, 2.05,
        2.1 , 2.15, 2.2 , 2.25, 2.3 , 2.35, 2.4 , 2.45, 2.5 , 2.55, 2.6 ,
        2.65, 2.7 , 2.75, 2.8 , 2.85, 2.9 , 2.95, 3.0, 3.05, 3.1 , 3.15, 3.2 , 3.25, 3.3 , 3.35, 3.4 , 3.45, 3.5])

#,
w_vec_ring = np.array([1.5,1.8,2.5])
meanratiovecvsw = np.zeros((len(w_vec_ring),3))
i = 0
#path = '/Users/maria/cluster_CSG2/Documents/project_phase_precession/LAING_setup'
os.chdir(path)
for omega in w_vec_ring:
    
    data1 = np.loadtxt('theta_A_0.95_N_100_beta_0.2_ic_chimera_om_'+str(omega)+'_tfinal_10000.txt')
    theta = data1[:,1:]
    deltat = 1000
    ndeltat = int(deltat/dt)
    mpv = np.zeros(N)
    
    length = len(theta)
    nl = int(length/ndeltat)
    meanratiovec = np.zeros((nl,2))
    for n in range(nl):
        for j in range(N):
            y = theta[n*ndeltat:(n+1)*ndeltat,j]
            peaks_t = find_peaks(y)
            peaks_1 = peaks_t[0]
            rotations_theta = len(peaks_1)
            peaks_1 = peaks_t[0]
            mean_phase_velocity_theta = 2*np.pi*(rotations_theta)/deltat
            mpv[j] = np.round(np.mean(mean_phase_velocity_theta),3)
        
        # RATIO
        sincro = np.min(mpv)
        mask_1 = (mpv == sincro)
        mask_2 = (mpv != sincro)
        pop_sincro = mpv[mask_1]
        pop_unsin = mpv[mask_2]
    
        ratio = sincro/pop_unsin
        
        meanratiovec[n,0] = np.mean(ratio)
        meanratiovec[n,1] = np.std(ratio) 
        
        if (n==1):
            np.savetxt('mpv_omega_'+str(omega)+'_N_'+str(N)+'.txt',mpv)

    meanratiovecvsw[i,0] = np.round(omega,2)
    meanratiovecvsw[i,1] = np.round(np.mean(meanratiovec[:,0]),3)
    meanratiovecvsw[i,2] = np.round(np.mean(meanratiovec[:,1]),3)
    i = i+1
    
np.savetxt('meanratio_vs_w_N_100_different_delta.txt',meanratiovecvsw)
    
