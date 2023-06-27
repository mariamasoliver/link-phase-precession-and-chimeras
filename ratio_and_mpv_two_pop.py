#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 18:29:36 2023

@author: maria
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy.signal import find_peaks



w_vec = np.array([1, 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.5 ,
       1.55, 1.6 , 1.65, 1.7 , 1.75, 1.8 , 1.85, 1.9 , 1.95, 2, 2.05,
       2.1 , 2.15, 2.2 , 2.25, 2.3 , 2.35, 2.4 , 2.45, 2.5 , 2.55, 2.6 ,
       2.65, 2.7 , 2.75, 2.8 , 2.85, 2.9 , 2.95, 3])

ratio = np.zeros(len(w_vec))


path = '/Users/maria/cluster_CSG/Documents/project_phase_precession/ic_0/'
os.chdir(path)
for i in range (len(w_vec)):
    print(w_vec[i])
    beta = 0.025
    K = 0.1
    TC = 0.012
    ic = 0
    nodes = 3 #change the number if you are working with N = 25
    llavor =3
    m = 2*nodes
    w = w_vec[i]
    
    """Removing transients and some extra"""
    t1 = 5000
    ti_file=t1
    dt = 0.1
    nti_file = int(ti_file/dt)

    data = np.loadtxt('theta_A_'+str(K)+'_N_'+str(nodes)+'_beta_'+str(beta)+'_ic_'+str(ic)+'_omega_'+str(w_vec[i])+'.txt',skiprows=nti_file)
    data2 = np.loadtxt('phi_A_'+str(K)+'_N_'+str(nodes)+'_beta_'+str(beta)+'_ic_'+str(ic)+'_omega_'+str(w_vec[i])+'.txt',skiprows=nti_file)
    
    dtfile = 0.1
    
    theta = data[:,1:]
    phi = data2[:,1:]
    
    ntf_output = len(theta[:,0])
    t_final_output = ntf_output*dtfile
    time = np.arange(0,t_final_output, dtfile)
    

    
    ndeltat = 10000


    deltat = ndeltat*dtfile
    ntotal_deltat = int(ntf_output/ndeltat)
    ndeltat_vec = np.arange(ntotal_deltat)*ndeltat
    
    """
    Rotations
    """
       
    j = 0
    rotations_theta = np.zeros(nodes)
    for l in range (nodes):
        peaks_t = find_peaks(theta[:ndeltat,l])
        peaks_1 = peaks_t[0]
        rotations_theta[l] = len(peaks_1)
    
    rotations_phi = np.zeros(nodes)
    for l in range (nodes):
        peaks_t = find_peaks(phi[:ndeltat,l])
        peaks_1 = peaks_t[0]
        rotations_phi[l] = len(peaks_1)
    
    """
    MPV theta i phi
    """
    
    mean_phase_velocity_theta = 2*np.pi*(rotations_theta)/deltat
    mean_phase_velocity_phi = 2*np.pi*(rotations_phi)/deltat
    
    
    mpvtheta = np.round(np.mean(mean_phase_velocity_theta),3)  #change the number if you are working with N = 25
    mpvphi = np.round(np.mean(mean_phase_velocity_phi),3)
    
    #the sincronized population has the smallest mean phase velocity
    #for the two pop chimera, the unsincro population has the same mpv for all oscillators, as
    #opposite to the chimera on a ring. 
    
    if (mpvphi[0] < mpvtheta[0]):
        mpvsin = mpvphi[0]
        mpvunsin = mpvtheta[0]
    else:
        mpvsin = mpvtheta[0]
        mpvunsin = mpvphi[0]
        

    mpv = np.append(mpvtheta,mpvphi)
    ratio = mpvsin/mpvunsin
    
    np.savetxt('mpv_A_'+str(K)+'_delta_'+str(deltat)+'_ic_0_w_'+str(w_vec[i])+'_TC_1_N_3.txt',mpv)
    
