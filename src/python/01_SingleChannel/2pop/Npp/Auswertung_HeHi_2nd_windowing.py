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

root='/home/raid3/tkunze/Documents/Programming/JansenModels/Kunze_JuR/003_ModelDevelopment/001_Unifying_Framework/EIPy_StateSpaceExploration/'

currpath = 'EI_Cluster/Data/2nd/2pop/2popNpp_HeHi_full2_windowing'
currpath = 'EI_Cluster/Data/2nd/2pop/2popNpp_widerImpulseRange'
os.chdir(root+currpath)

###############################################################################
#
# Function definitions
#
###############################################################################

def cleargrid2(state_grid):
    [x,y,z]=np.shape(state_grid)
    out_grid=np.zeros((x,y))    
    for i in range(x):
        for j in range(y):
            if state_grid[i,j,0] > 0.004:           # (1,min) Region was above
            
                if state_grid[i,j,4] > 0.004:       # (3,min) Region is above
                
                    if state_grid[i,j,2] > 0.004:       # (2,min) region was temporary above
                        out_grid[i,j] = 11                  # was above, (stayed up) and is now above

                    elif state_grid[i,j,2] < 0.004:     # (2,min) region was temporary below
                        out_grid[i,j] = 10                  # was above, was temporary below and is now above again
                        
                elif (state_grid[i,j,5] > 0.004) and (state_grid[i,j,4] < 0.004):   # (3,min,max) Region is oscillating
                    out_grid[i,j] = 9                         # was above, is now oscillating

                elif state_grid[i,j,4] < 0.004:     # (3,min) Region is now low
                    out_grid[i,j] = 8                       # was above, became below
                    
            elif (state_grid[i,j,1] > 0.004) and (state_grid[i,j,0] < 0.004):   # (1,min,max) Region was oscillating       
                if state_grid[i,j,4] > 0.004:       # (3,min) Region is above
                    out_grid[i,j] = 7                       # was oscillating, is now above
                
                elif (state_grid[i,j,5] > 0.004) and (state_grid[i,j,4] < 0.004):   # (3,min,max) Region is oscillating
                    out_grid[i,j] = 6                         # was oscillating, is now oscillating
                
                elif state_grid[i,j,4] < 0.004:       # (3,min) Region is below
                    out_grid[i,j] = 5                         # was oscillating, is now below
            
            elif state_grid[i,j,1] < 0.004:         # (1,max) Region was below
            
                if state_grid[i,j,4] > 0.004:       # (3,min) Region is above
                    out_grid[i,j] = 4                   # was below, and is now above
                    
                elif (state_grid[i,j,5] > 0.004) and (state_grid[i,j,4] < 0.004):   # (3,min,max) Region is oscillating
                    out_grid[i,j] = 3                   # was below, and is now oscillating
                    
                elif state_grid[i,j,5] < 0.004:     # (3,max) Region is now below
                    
                    if state_grid[i,j,3] > 0.004:   # (2,max) Region was temporary above
                        out_grid[i,j] = 2               # was below, was temporary above and is now below
                    
                    elif state_grid[i,j,3] < 0.004: # (2,max) Region wasn't temporary up
                        out_grid[i,j] = 1               # was below and stayed permantly below
                    
            else:
                raise ValueError('2')
                            
    return out_grid       

###############################################################################
#
# Analysis
#
###############################################################################

# baselevel, 2: endlevel, 3:max()
from matplotlib import cm   
import matplotlib.pyplot as plt


#full2
hirange = np.arange(2,8.1,0.5)*1e-3
herange = np.arange(7.,0.99,-0.5)*1e-3
length_range = np.arange(500,1501,10)
intensity_range = np.arange(250,49,-2)

#wide Impulse range
hirange = np.arange(2,8.1,0.5)*1e-3
herange = np.arange(7.,0.99,-0.5)*1e-3
length_range = np.arange(10,1501,10)
intensity_range = np.arange(250,-1,-2)


low_lenge=np.min(length_range)
high_lenge=np.max(length_range)
low_inte=np.min(intensity_range)
high_inte=np.max(intensity_range)




#%%

state_grid=np.zeros((len(herange),len(hirange),len(intensity_range),len(length_range),6))     
i=0
j=0
for he in herange:
    j=0    
    for hi in hirange:
        state_grid[i,j,:,:,:]=np.load('RUM_Dec_meas_Npp113k4_HeHi_le500t1500i10msInt50t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000))
        #state_grid[i,j,:,:,:]=np.load('RUM_Dec_meas_Npp113k4_le10t1500i10msInt0t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000))        
        print "he:%.2f | hi:%.2f " %(he*1000,hi*1000)
        j=j+1
    i=i+1


glob_low_val_base=np.min(state_grid[:,:,:,:,0])
glob_high_val_base=np.max(state_grid[:,:,:,:,1])

glob_low_val_max=np.min(state_grid[:,:,:,:,2])
glob_high_val_max=np.max(state_grid[:,:,:,:,3])

glob_low_val_dete=np.min(state_grid[:,:,:,:,4])
glob_high_val_dete=np.max(state_grid[:,:,:,:,5])



scale_base=np.zeros((2,2))
scale_base[0]=glob_low_val_base
scale_base[1]=glob_high_val_base
scale_detec=np.zeros((2,2))
scale_detec[0]=glob_low_val_dete
scale_detec[1]=glob_high_val_dete
scale_max=np.zeros((2,2))
scale_max[0]=glob_low_val_max
scale_max[1]=glob_high_val_max
scale_colorcode=np.array([[11,11],[10,10],[9,9],[8,8],[7,7],[6,6],[5,5],[4,4],[3,3],[2,2],[1,1]])


fig=plt.figure(100)
plt.clf()
plt.subplot(141)
plt.imshow(scale_base, aspect='auto', interpolation='none')
plt.colorbar()
plt.title('scale base level')
plt.subplot(142)
plt.imshow(scale_detec, aspect='auto', interpolation='none')
plt.colorbar()
plt.title('scale detector level')
plt.subplot(143)
plt.imshow(scale_max, aspect='auto', interpolation='none')
plt.colorbar()
plt.title('scale max level')
plt.subplot(144)
plt.imshow(scale_colorcode, aspect='auto', interpolation='none',cmap=cm.Accent)
plt.colorbar()
plt.title('scale color code')


#%%

plt.figure(3) 
plt.clf()
plt.figure(4) 
plt.clf()          
i=0
j=0
n=1

for he in herange:
    j=0    
    for hi in hirange:
        

        print "he:%.2f | hi:%.2f " %(he*1000,hi*1000)
        
        #if np.min(state_grid[:,:,1]) < glob_low_val: 
        #    glob_low_val=np.min(state_grid[:,:,1])
        #    print he,hi,glob_low_val
        #if np.max(state_grid[:,:,1]) > glob_high_val: 
        #    glob_high_val=np.max(state_grid[:,:,1])
        #    print he,hi,glob_high_val,1
        plt.figure(3)
        plt.subplot(len(herange),len(hirange),n)
        plt.text(0.5, 0.5, "'he:%.2f,hi:%.2f,#:%.0f" %(he*1000,hi*1000,n), size=8, rotation=0.,
        ha="center", va="center",
        bbox=dict(boxstyle="square",
                  ec=(1., 1., 1.),
                  fc=(1., 1., 1.),
                  )
         )
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])

    
        plt.figure(4)
        plt.subplot(len(herange),len(hirange),n)
        #current_grid=state_grid[i,j,:,:,:]
        current_grid=cleargrid2(state_grid[i,j,:,:,:])
        #sectionedgrid=classifygrid(state_grid)
        current_grid[-1,-1]=11
        current_grid[1,2]=1
        plt.imshow(np.flipud(current_grid[:,:]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none',cmap=cm.Accent)
        #plt.title('he:%.1f,hi:%.0f' %(he*1000,hi*1000))
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])
        if i==0:
            plt.title('hi:%.2f' %(hi*1000))
        if np.mod(i,2)==0 and j==0:
            plt.ylabel('%.2f' %(he*1000))
        j=j+1
        n=n+1
        if n==3000:
            raise ValueError('dd')

    i=i+1

plt.figure(2)
plt.tight_layout()
plt.figure(4)
plt.tight_layout()
#cb=plt.colorbar()
#plt.ylabel('intensity')
#plt.xlabel('length')
print "fertig"        

def isosci(state_grid):
    [x,y,z]=np.shape(state_grid)
    out_grid=np.zeros((x,y)) 
    for i in range(x):
        for j in range(y):
            if state_grid[i,j,1]-state_grid[i,j,0] > 0.000000000001:           # region is oscillating
                out_grid[i,j] = 1
                
    return out_grid
    
n=1
j=0
i=0
plt.figure(6)
plt.clf()
current_grid=np.zeros((len(intensity_range),len(length_range)))
for he in herange:
    j=0    
    for hi in hirange:
        plt.figure(6)
        plt.subplot(len(herange),len(hirange),n)
        #current_grid=isosci(state_grid[i,j,:,:,:])
        #current_grid[-1,-1]=0
        #current_grid[1,2]=1
        currentgrid=state_grid[i,j,:,:,1]-state_grid[i,j,:,:,0]
        print current_grid[12,12]
        current_grid[-1,-1]=0
        #current_grid[-1,-2]=glob_high_val_base
        #plt.imshow(np.flipud(current_grid[:,:]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none',cmap=cm.Accent)
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])
        if i==0:
            plt.title('hi:%.2f' %(hi*1000))
        if np.mod(i,2)==0 and j==0:
            plt.ylabel('%.2f' %(he*1000))
        j=j+1
        n=n+1
    i=i+1

plt.figure(6)
plt.tight_layout()    
#%%

"""
single fingerprint plots
"""
he=2.5e-3
hi=5e-3
a=np.load('RUM_Dec_meas_Npp113k4_HeHi_le500t1500i10msInt50t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000))

plt.figure(2)
plt.clf()
current_grid=cleargrid2(a[:,:,:])
current_grid[-1,-1]=11
current_grid[0,0]=1
plt.imshow(np.flipud(current_grid[:,:]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none',cmap=cm.Accent)

currpath=os.getcwd()
os.chdir('/home/raid3/tkunze/Documents/kurzzeitig')
plt.savefig('fingerprint_2popNpp_he%.2f_hi%.0f.eps' %(he*1000,hi*1000), format='eps', dpi=1000)
os.chdir(currpath)

#%%
he=7e-3
hi=3.5e-3
a=np.load('RUM_Dec_meas_Npp100_HeHi_le500t1500i10msInt50t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000))

base=np.flipud(a[:,:,0])
base[-1,-1]=glob_low_val_base
base[-1,-2]=glob_high_val_base

detector=np.flipud(a[:,:,1])
detector[-1,-1]=glob_low_val_dete
detector[-1,-2]=glob_high_val_dete

maxlevel=np.flipud(a[:,:,2])
maxlevel[-1,-1]=glob_low_val_max
maxlevel[-1,-2]=glob_high_val_max

plt.figure(14)
plt.clf()
plt.imshow(base, aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none')
plt.colorbar()            
plt.title('baselevel, he%.1f:, hi:%.1f' %(he*1000,hi*1000))

plt.figure(15)
plt.clf()
plt.imshow(detector, aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none')
plt.colorbar()            
plt.title('detektorlevel, he%.1f:, hi:%.1f' %(he*1000,hi*1000))

plt.figure(16)
plt.clf()
plt.imshow(maxlevel, aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none')
plt.colorbar()            
plt.title('maxlevel, he%.1f:, hi:%.1f' %(he*1000,hi*1000))
            



            
#%%
currpath = '/home/raid3/tkunze/Documents/Programming/JansenModels/Kunze_JuR/003_ModelDevelopment/001_Unifying_Framework/EIPy_StateSpaceExploration/EI_Cluster/Figures'
os.chdir(currpath)

plt.figure(2)
plt.savefig('RUM_Detektor_2popNpp_sections.pdf', format='pdf', dpi=1000)
plt.figure(100)
plt.savefig('RUM_Detektor_2popNpp_colorcode.pdf', format='pdf', dpi=1000)
plt.figure(3)
plt.savefig('RUM_Detektor_2popNpp_1_detectorlevel.pdf', format='pdf', dpi=1000)
plt.figure(4)
plt.savefig('RUM_Detektor_2popNpp_2_wideImpulseRange.pdf', format='pdf', dpi=1000)
plt.figure(5)
plt.savefig('RUM_Detektor_2popNpp_1_sectioned.pdf', format='pdf', dpi=1000)
plt.figure(14)
plt.savefig('RUM_Detektor_2popNpp_1_baselevel_single.pdf', format='pdf', dpi=1000)
plt.figure(15)
plt.savefig('RUM_Detektor_2popNpp_1_detelevel_single.pdf', format='pdf', dpi=1000)
plt.figure(16)
plt.savefig('RUM_Detektor_2popNpp_1_maxlevel_single.pdf', format='pdf', dpi=1000)
#plt.close('all')
