# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 17:15:03 2015

@author: Tim Kunze

Copyright (C) 2015, Tim Kunze. All rights reserved.



This script is to check whether there are fix points or limit cycles in the non activated state (pext=0s-1)
with in the respective range of He-Hi-diagram (3pop)

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
import matplotlib.pyplot as plt

import sys
sys.path.append("/home/raid3/tkunze/Documents/Programming/JansenModels/Kunze_JuR/003_ModelDevelopment/001_Unifying_Framework") 


#2pop version
import Models.Generic_fin_01 as FCV2

import Simulation_And_Analysis.Sim_Simulation_003 as simulate
  

currpath = '/home/raid3/tkunze/Documents/Programming/JansenModels/Kunze_JuR/003_ModelDevelopment/001_Unifying_Framework/EIPy_StateSpaceExploration'
os.chdir(currpath)


#%% 
###############################################################################
#
# Main
#
###############################################################################


JR = FCV2.JuR()
dt = 1000e-6/0.01       # to account for implementation of the model in dimensionless form
JR.integ_stepsize = dt
JR.b1=0		#controls connection pe and ep: 1-> connected
JR.b2=0		#controls input               : 1-> input to EI
JR.b3=0		#controls self-conn of PC     : 1-> No selfconn PP
JR.b4=0		#controls self-conn of IIN    : 1-> No selfconn II
JR.n=2
JR.coupling_II = np.zeros((2,2))
JR.coupling_EX = np.zeros((2,2))
JR.distanceMatrix = (np.ones((2,2)) - np.identity(2))*0.001
JR.init = np.zeros((10,JR.n))           
JR.c_In_ii=0
JR.c_In_ex=0
JR.c_pp=113.4
JR.c_ii=33.25
JR.configure()

#%%

###############################################################################
#
## RUM show case
###############################################################################

#%% 2pop version
lenge=940
inte=100

JR.H_e=7e-3
JR.H_i=15e-3 

t_simulation = 5/0.01
N=t_simulation/dt
time = np.arange(0,N*dt,dt)
 

p_sim_ex = np.zeros((N,JR.n))
p_sim_i  = np.zeros((N,JR.n))
p_sim_ex[1000:1000+lenge,:]=inte

signal,sig_ei,sig_ii,impact,data = simulate.simulate_network_SHC(JR,p_sim_ex,p_sim_i,t_simulation)

plt.figure(25)
plt.clf()
plt.subplot(211)
plt.plot(time*0.01,signal[:,0]/560)
plt.grid()
plt.plot([(1000+lenge)/1000, (1000+lenge)/1000 ],[-0.01, 0.015],'k')
plt.xlabel('time in s')
plt.ylabel('pyrapot in V')
plt.title('he%.1f:, hi:%.1f | inte=%.1f | length:%.1fms' %(JR.H_e*1000,JR.H_i*1000,inte,lenge))
plt.subplot(212)
plt.plot(time*0.01,p_sim_ex+p_sim_i)
plt.ylabel('pext')
plt.xlabel('time in s')


###############################################################################
#
## IOsci test
###############################################################################

#%% 3pop version


t_simulation = 1.5/0.01

N=t_simulation/dt
time = np.arange(0,N*dt,dt)




p_sim_ex = np.zeros((N,JR.n))
p_sim_i  = np.zeros((N,JR.n))

he_range = np.arange(0,16.1,0.1)*1e-3
hi_range = np.arange(0,35.1,0.1)*1e-3

state_grid = np.zeros((len(he_range),len(hi_range),2))
sim_grid=np.zeros((len(he_range),len(hi_range),N))
i=0
j=0
for he in he_range:
    j=0
    JR.H_e=he
    for hi in hi_range:
        JR.H_i=hi
        print "he: %.2fmV | hi: %.2fmV" %(he*1000,hi*1000)

        signal,sig_ei,sig_ii,impact,data = simulate.simulate_network_SHC(JR,p_sim_ex,p_sim_i,t_simulation)
        state_grid[i,j,0] = np.min(signal[1000:,0])
        state_grid[i,j,1] = np.max(signal[1000:,0])
        sim_grid[i,j,:]=signal[:,0]
        j+=1
        
    i+=1

np.save('RUM_oscitest_2popNii_he0i0k1t16_hi0i0k1t35_pext0_state_grid.npy',state_grid)
np.save('RUM_oscitest_2popNii_he0i0k1t16_hi0i0k1t35_pext0_sim_grid.npy',sim_grid)




#
###############################################################################
#
## Analysis
#
###############################################################################

stategrid=np.load('RUM_oscitest_2popNpp_he0i0k1t16_hi0i0k1t20_pext0_state_grid.npy')

simgrid=np.load('RUM_oscitest_2popNpp_he0i0k1t16_hi0i0k1t20_pext0_sim_grid.npy')
he_pos=80
hi_pos=100
plt.figure(50)
plt.clf()
plt.plot(simgrid[he_pos,hi_pos,:]/560)


[rows,cols,vals]=np.shape(stategrid)

oscigrid=np.zeros((rows,cols))
for i in range(rows):
    for j in range(cols):
        if (stategrid[i,j,1]/560-stategrid[i,j,0]/560)<0.001:
            oscigrid[i,j]=0                         #no oscillations
        else:
            oscigrid[i,j]=1                         # oscillations occur
 

he_range = np.arange(0,16.1,0.1)*1e-3
hi_range = np.arange(0,35.1,0.1)*1e-3

min_hi=np.min(hi_range)
max_hi=np.max(hi_range)
min_he=np.min(he_range)
max_he=np.max(he_range)

plt.figure(2)
#state_grid=cleargrid(state_grid)
plt.clf()
plt.imshow(np.flipud(oscigrid), aspect='auto', extent = (min_hi,max_hi,min_he,max_he),interpolation='none')
plt.ylabel('He in V')
plt.xlabel('Hi in V')
plt.title('Oscillation diagram')
plt.colorbar()

plt.savefig('RUM_2popNii_oscillationstest.pdf', format='pdf', dpi=1000)
#




