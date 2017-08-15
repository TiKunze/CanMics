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

Npp ändern!!!!!

import numpy as np
import sys
import scipy as sc
import os               # to enable some C commands (cwd,listdir)


currpath = '/usr/wrk/people9/tiku2449/EI_RUM/001_Unifying_Framework/RUM_Exploration/2pop'
os.chdir(currpath)

import sys
sys.path.append("/usr/wrk/people9/tiku2449/EI_RUM/001_Unifying_Framework") 


import Models.Generic_fin_01 as FCV
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


dt = 1000e-6/0.01       # to account for implementation of the model in dimensionless form

JR = FCV.JuR()

JR.integ_stepsize = dt

JR.b1=0		#controls connection pe and ep: 1-> connected
JR.b2=0		#controls input               : 1-> input to EI
JR.b3=0		#controls self-conn of PC     : 1-> No selfconn PP
JR.b4=1		#controls self-conn of IIN    : 1-> No selfconn II

JR.n=2

JR.coupling_II = np.zeros((2,2))
JR.coupling_EX = np.zeros((2,2))
JR.distanceMatrix = (np.ones((2,2)) - np.identity(2))*0.001


JR.init = np.zeros((10,JR.n))           
JR.c_In_ii=0
JR.c_In_ex=0
JR.configure()


#%%
###############################################################################
#
## Activation Diagram RUM with modulation of input to II
#
###############################################################################

t_simulation = 5/0.01

N=t_simulation/dt
time = np.arange(0,N*dt,dt)

JR.H_e=he
JR.H_i=hi

length_range = np.arange(500,1501,10)
intensity_range = np.arange(50,251,2)
state_grid = np.zeros((len(intensity_range),len(length_range),6))  



p_sim_ex = np.zeros((N,JR.n))
p_sim_i  = np.zeros((N,JR.n))

i=0
j=0
for ins in intensity_range:
    j=0
    for le in length_range:
        p_sim_ex = np.zeros((N,JR.n))      
        p_sim_ex[1000:1000+le,:]=ins
        signal,sig_ei,sig_ii,impact,data = simulate.simulate_network_SHC(JR,p_sim_ex,p_sim_i,t_simulation)
        state_grid[i,j,0] = np.min(signal[500:999,0]/560)
        state_grid[i,j,1] = np.max(signal[500:999,0]/560)         
        state_grid[i,j,2] = np.min(signal[1100:3500,0]/560)
        state_grid[i,j,3] = np.max(signal[1100:3500,0]/560)      
        state_grid[i,j,4] = np.min(signal[4000:,0]/560)
        state_grid[i,j,5] = np.max(signal[4000:,0]/560)
        print "len: %.0f | int: %.0f | he: %.2fmV | hi: %2.fmV" %(le, ins, he*1000,hi*1000)
        j+=1
        
    i+=1


np.save('RUM_Dec_meas_Npp113k4_HeHi_le500t1500i10msInt50t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000),state_grid)
