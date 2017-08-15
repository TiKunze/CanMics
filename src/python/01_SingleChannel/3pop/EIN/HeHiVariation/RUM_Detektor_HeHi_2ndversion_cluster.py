# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 17:15:03 2015

@author: Tim Kunze

Copyright (C) 2015, Tim Kunze. All rights reserved.



This script is a modified version of the RUM Detector:

instead of sweeping over He and Hi in every diagram, we sweep over lenge and intensity of the impulse (as in the actiation plot)

"""

###############################################################################
#
# Imports
#
###############################################################################



import numpy as np
import sys
import scipy as sc
import os               # to enable some C commands (cwd,listdir)


currpath = '/usr/wrk/people9/tiku2449/EI_RUM/001_Unifying_Framework/RUM_Exploration'
os.chdir(currpath)

import sys
sys.path.append("/usr/wrk/people9/tiku2449/EI_RUM/001_Unifying_Framework") 


import Models.JuRClass_fin_006 as FCV
import Simulation_And_Analysis.Sim_Simulation_003 as simulate


while len(sys.argv) > 1:
    option = sys.argv[1];                                             del sys.argv[1]
    if   option == '-he':  he = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-hi':  hi = float(sys.argv[1].replace(',','.'));  del sys.argv[1]   
    else:
        print 'Options invalides :',option,'->',sys.argv[0]


#%% 
###############################################################################
#
# Main
#
###############################################################################


dt = 1000e-6

JR = FCV.JuR()

JR.integ_stepsize = dt

JR.n=2
JR.coupling = np.array([[0.0,0.0],[0.0,0.0]])   # 

JR.distanceMatrix = np.array([[0.0,0.01],[0.0,0.0]]) # important!!
JR.init = np.zeros((8,JR.n))           
JR.c_e=0        # only relevant for connected areas
JR.c_i=0        # only relevant for connected areas
JR.c_py=30       # only relevant for connected areas
JR.configure()


#%%
###############################################################################
#
## Activation Diagram RUM with modulation of input to II
#
###############################################################################

t_simulation = 5

N=t_simulation/dt
time = np.arange(0,N*dt,dt)

JR.H_e=he
JR.H_i=hi

p_sim_py  = np.zeros((N,JR.n))
p_sim_e  = np.zeros((N,JR.n))
p_sim_i  = np.zeros((N,JR.n))
length_range = np.arange(500,1501,10)
intensity_range = np.arange(50,251,2)
state_grid = np.zeros((len(intensity_range),len(length_range),3))    

i=0
j=0
for ins in intensity_range:
    j=0
    for le in length_range:
        p_sim_e  = np.zeros((N,JR.n))        
        p_sim_e[1000:1000+le,:] = ins
        signal,sig_ei,sig_ii,impact,data = simulate.simulate_network_006(JR,p_sim_py,p_sim_e,p_sim_i,t_simulation)
        state_grid[i,j,0] = np.mean(signal[999,0])        
        state_grid[i,j,1] = np.mean(signal[4000:,0])
        state_grid[i,j,2] = np.max(signal[900:,0])
        print "len: %.0f | int: %.0f | he: %.2fmV | hi: %2.fmV" %(le, ins, he*1000,hi*1000)
        j+=1
        
    i+=1


#dataa=length_range,intensity_range,state_grid
np.save('RUM_Dec_meas_full2_le500t1500i10msInt50t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000),state_grid)
#np.save('RUM_Dec_sim_le500t1500i10msInt70t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000),signal)
#np.save('RUM_Dec_data_le500t1500i10msInt70t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000),dataa)



#
#
#def cleargrid(state_grid):
#    [x,y,z]=np.shape(state_grid)
#    for i in range(x):
#        for j in range(y):
#            if state_grid[i,j,1] > 0.004:
#                state_grid[i,j,1] = 0.006
#            elif state_grid[i,j,1] < 0.004:
#                state_grid[i,j,1] = -0.002
#            else:
#                raise ValueError('Error')
#                print "ERROR"
#                
#    return state_grid
###            
#%% Analysis
#import matplotlib.pyplot as plt
#hirange = np.arange(19,26,1)*1e-3
#herange = np.arange(2.5,4.1,0.25)*1e-3
#    
#
#glob_low_val=1e3
#glob_high_val=-1e3
#
#for he in herange:
#    for hi in hirange:
#        a=np.load('RUM_Detector2_Imple500t1500i10msInt70t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000))
#        
#          
#        length_range=a[0]
#        intensity_range=a[1]
#        state_grid=a[2]
#        
#        low_lenge=np.min(length_range)
#        high_lenge=np.max(length_range)
#        low_inte=np.min(intensity_range)
#        high_inte=np.max(intensity_range)
#        
#        if np.min(state_grid[:,:,1]) < glob_low_val: 
#            glob_low_val=np.min(state_grid[:,:,1])
#            print he,hi,glob_low_val
#        if np.max(state_grid[:,:,1]) > glob_high_val: 
#            glob_high_val=np.max(state_grid[:,:,1])
#            print he,hi,glob_high_val,1
#        
#        plt.figure(2) 
#        plt.clf()
#        state_grid=cleargrid(state_grid)
#        plt.imshow(np.flipud(state_grid[:,:,1]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none')
#        plt.ylabel('intensity')
#        plt.xlabel('length')
#        plt.title('Detektor Diagram,he:%.0fms, hi:%.0fpps' %(he*1000,hi*1000))
#        cb=plt.colorbar()
#        plt.savefig('RUM_Detektor2_Imple500t1500i10msInt70t250i2_He%.2fmV_Hi%.1fmV.pdf'  %(he*1000,hi*1000), format='pdf', dpi=1000)
#        plt.close()
#        #
#        
#        # baselevel plot hier zwecklos, da baselevel bei allen stimuli gleich
#        plt.figure(2) 
#        plt.clf()
#        #state_grid=cleargrid(state_grid)
#        plt.clf()
#        plt.imshow(np.flipud(state_grid[:,:,0]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none')
#        plt.ylabel('intensity')
#        plt.xlabel('length')
#        plt.title('Detektor Diagram,Baselevels,he:%.0fmV, hi:%.0fmV' %(he*1000,hi*1000))
#        plt.colorbar()
#        #plt.savefig('RUM_Detektor_Baselevel_Imple%.0fmsInt%.0f_He2.5t7.0i0k05_Hi10t25i0k1_1.pdf' %(lenge,inte), format='pdf', dpi=1000)
#        
#        plt.close('all')
