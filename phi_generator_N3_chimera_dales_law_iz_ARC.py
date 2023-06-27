#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recurrent Neural Network of Spiking Neurons
Force Method with Izhikevich Network 
"""
import numpy as np 
def main(path, N, q):    
    
    import time
    it = time.time()
    
    T = 100000 #Total time in ms
    dt = 0.04 #Integration time step in ms 
    nt = int(T/dt) #Time steps

    # Izhikevich Parameters
    C = 250  #capacitance
    vr = -60   #resting membrane 
    b = -2  #resonance parameter 
    ff = 2.5  #k parameter for Izhikevich, gain on v 
    vpeak = 30  # peak voltage
    vreset = -65 # reset voltage 
    vt = vr+40-(b/ff) #threshold 
    u = np.zeros(N)  #initialize adaptation 
    a = 0.01 #adaptation reciprocal time constant 
    d = 200 #adaptation jump current 
    tr = 2  #synaptic rise time 
    td = 20 #decay time 

    G =5*10**3 #Gain on the static matrix with 1/sqrt(N) scaling weights.  
    Q =int(0.1*q*10**3) #Gain on the rank-k perturbation modified by RLS.
    
    
    #
    mu = 0.5 #ratio of E to total 
    NE = int(mu*N)

    
    p = 0.1 #network sparsity
    p_temp = p/mu
    
    #Storage variables for synapse integration  
    IPSC = np.zeros(N) #post synaptic current 
    h = np.zeros(N)
    r = np.zeros(N)
    hr = np.zeros(N)
    JD = np.zeros(N)
    
    #-----Initialization---------------------------------------------
    v = vr+(vpeak-vr)*np.random.rand(N) #initial distribution 
    v_ = v
    
    
    #############
    #supervisor
    nodes = 3
    m = 2*2*nodes
    K = 0.1
    beta = 0.025
    TC = 0.012
    ic = 0
    #Removing transients
    ti_file = 2000
    nti_file = int(ti_file/dt)
    data = np.loadtxt('theta_A_'+str(K)+'_N_'+str(nodes)+'_beta_'+str(beta)+'_ic_'+str(ic)+'_TC_'+str(TC)+'_omega_2.57.txt',skiprows=nti_file)
    data2 = np.loadtxt('phi_A_'+str(K)+'_N_'+str(nodes)+'_beta_'+str(beta)+'_ic_'+str(ic)+'_TC_'+str(TC)+'_omega_2.57.txt',skiprows=nti_file)
       
    
    theta = np.transpose(data[nti_file:,1:])
    phi = np.transpose(data2[nti_file:,1:])
    
    ntf_file = len(theta[0,:])
    
    sup = np.zeros((m,ntf_file))
    sup[0:nodes] = np.sin(theta)
    sup[nodes:2*nodes] = np.cos(theta)
    sup[2*nodes:3*nodes] = np.sin(phi)
    sup[3*nodes:4*nodes] = np.cos(phi)
    
    OMEGA = np.zeros((N,N))
    np.random.seed(3)
    for i in range(0,N):
        A = np.random.normal(0., 1.0, (N))
        OMEGA[i,0:NE]=np.multiply(A[0:NE],A[0:NE]>0)
        OMEGA[i,NE:N]=np.multiply(A[NE:N],A[NE:N]<0)
    
    
    B = np.random.rand(N,N)<p_temp #NxN Matrix with a uniform distribution. 
                              #It will return true and false, depending on p. 
    OMEGA = G*np.multiply(OMEGA,B)/np.sqrt(N*p) #Initial weight matrix
    
    #Kill EE and EI connections.  
    OMEGA[:N,:NE] = 0
    
    
    eta = (2*np.random.rand(N,m)-1)*Q #random theta variable ranging from -1,1
    etaminus = np.zeros((N,m))
    etaplus = eta + etaminus
    etaplus[etaplus<0] = 0
    etaminus = eta - etaplus
    
    output = np.zeros(m)  #initial approximant
    output1 = np.zeros(m)
    output2 = np.zeros(m)
    
    #initialize decoders & RLS method
    dec = np.zeros((N,m)) #initial decoder.  Best to keep it at 0.  
    dec1 = np.zeros((N,m))
    dec2 = np.zeros((N,m))
    
    tmax = T #max saving time for the neurons that fired
    ntmax = int(tmax/dt)
     
    n_fired = np.zeros((ntmax,N))  #If you want to store spike times, 
    
    BIAS = np.zeros(N) 
    BIAS[:NE]= -20;  #Background current to Excitatory neurons, to shut them off
    BIAS[NE:]= 3000; #Background current to Inhibitory neurons 
    
    
    
    ## 
    la = 0.5
    Pinv = np.eye(N)/la #initial correlation matrix, coefficient is the regularization constant as well 
    step_rls = 20 #optimize with RLS only every 50 steps 
    step_plot = 2 #not saving every time step
    imin = int(0.1*T/dt) #time before starting RLS, gets the network to chaotic attractor 
    icrit = int(0.8*T/dt) #end simulation at this time step 
    ntplot = int(nt/step_plot)
    
    save_output = np.zeros((ntplot,m))  #store the approximant 
    save_sup = np.zeros((ntplot,m))  #store the supervisor as used here
    save_dec = np.zeros((ntplot,5,m)) #store the decoders 
    
    save_v_exc = np.zeros((ntplot,10)) #Store voltage and adaptation variables for plotting 
    save_v_inh = np.zeros((ntplot,10)) #Store voltage and adaptation variables for plotting 
    
    ## SIMULATION 
    iplot = 0
    for i in range(nt): 
    ## EULER INTEGRATE
        I = IPSC + np.dot(etaplus,output1) + np.dot(etaminus,output2)  + BIAS  #postsynaptic current 
        v = v + dt*((ff*np.multiply((v-vr),(v-vt)) - u + I))/C  # v(t) = v(t-1)+dt*v'(t-1)
        u = u + dt*(a*(b*(v_-vr)-u)) #same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
        ## 
        index = (v>=vpeak)
        if (len(v[index])>0):
            JD = OMEGA[:,index].sum(axis=1) #compute the increase in current due to spiking  
            if(i<ntmax):
                n_fired[i,:] = index  #comment this if you want to store spike times.  
        
        #synapse for double np.exponential        
        IPSC = IPSC*np.exp(-dt/tr) + h*dt
        h = h*np.exp(-dt/td) + JD*(len(v[index])>0)*1/(tr*td)  #Integrate the current        
        r = r*np.exp(-dt/tr) + hr*dt 
        hr = hr*np.exp(-dt/td) + (v>=vpeak)/(tr*td)
    
        output = np.dot(np.transpose(dec),r)                
        ##RLS
        if (i > imin):
            if (i < icrit):
                if (i%step_rls == 0):
                    e = output - sup[:,i]
                    q = np.dot(Pinv,r)               
                    Pinv = Pinv - (np.outer(q,np.transpose(q)))/(1+np.dot(np.transpose(r),q))
                    dec = dec - np.outer(np.dot(Pinv,r),e)
                    
                    decep = np.multiply(dec[0:NE,:],dec[0:NE,:]>0)
                    decem = np.multiply(dec[0:NE,:],dec[0:NE,:]<0)
                    decip = np.multiply(dec[NE:N,:],dec[NE:N,:]>0)
                    decim = np.multiply(dec[NE:N,:],dec[NE:N,:]<0)
                    
                    dec1[0:NE] = decep
                    dec1[NE:N] = decim
                    dec2[0:NE] = decem*mu/(1-mu)
                    dec2[NE:N] = decip*mu/(1-mu)
        
        u = u + d*(v>=vpeak)  #implements set u to u+d if v>vpeak, component by component. 
        v = v+np.multiply((vreset-v),(v>=vpeak)) #implements v = vreset if v>vpeak add 0 if false, add vreset-v if true, v+vreset-v = vreset
        v_ = v  # sets v(t-1) = v for the next itteration of loop
    
        ## Store  
        if (i%step_plot == 0):
            save_v_inh[iplot,:] = v[-10:]  
            save_v_exc[iplot,:] = v[:10]  
            save_output[iplot,:] = np.transpose(output)
            save_dec[iplot, :, :] = dec[:5,:]
            save_sup[iplot,:] = sup[:,i]
            iplot = iplot +1
    
    n_fired = n_fired.astype(int)
    n_fired_tf = n_fired == 1
    n_time_vec = np.arange(0,ntmax)
    # get the exact firing time
    n_fired_times = []
    for n in range(N):
        ar_temp = n_time_vec[n_fired_tf[:,n]]
        n_fired_times.append(ar_temp) #each entry is a list of the firing times for that neuron
    
    np.savetxt('v_chimera_Q_'+str(Q)+'_N_'+str(N)+'_T_'+str(T)+'_dl.txt', save_v_inh)
    np.savetxt('u_chimera_Q_'+str(Q)+'_N_'+str(N)+'_T_'+str(T)+'_dl.txt', save_v_exc)
    np.savetxt('chimera_Q_'+str(Q)+'_N_'+str(N)+'_T_'+str(T)+'_dl.txt', save_output)
    np.savetxt('sup_chimera_'+str(Q)+'_N_'+str(N)+'_T_'+str(T)+'_dl.txt',save_sup)
    
    np.save('tspike'+str(Q)+'_N_'+str(N)+'_T_'+str(T)+'_dl.npy',n_fired_times)
    np.savetxt('dec_'+str(Q)+'_N_'+str(N)+'_T_'+str(T)+'.txt', dec)
    
    ft = time.time()   
    time_elapsed = ft-it
    te = np.array([time_elapsed])
    np.savetxt('elapsed_time_Q_'+str(Q)+'_N_'+str(N)+'_T_'+str(T)+'.txt',te)




if __name__ == '__main__':
    import sys, os

    path = sys.argv[1]
    N = int(sys.argv[2])
    q = int(sys.argv[3])
    if not os.path.exists(path):
        print('Directory:', path, '\ndoes not exist, creating new one.')
        os.makedirs(path) 
    else:
        print('Saving to', path)
    main(path, N, q)

    
    
    
    
    
    
    
    
    
