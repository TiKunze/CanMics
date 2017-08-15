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

currpath = 'EI_Cluster/Data/2nd/2pop/2popNpp_HeHi_full1'

root='/home/raid3/tkunze/Documents/Programming/JansenModels/Kunze_JuR/003_ModelDevelopment/001_Unifying_Framework/EIPy_StateSpaceExploration/'
currpath = 'EI_Cluster/Data/2nd/3pop/full_2'

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
            if state_grid[i,j,0] > 0.004:           # (1) Region was up
            
                if state_grid[i,j,1] > 0.004:       # (3) Region is up
                    out_grid[i,j] = 5                   # was up, (stayed up) and is now up
                    
                elif state_grid[i,j,1] < 0.004:     # (3) Region is now low
                    
                    if state_grid[i,j,2] > 0.004:   # (2) Region was temporary up
                        out_grid[i,j] = 4               # was up, was temporary up and is now low
                    
                    elif state_grid[i,j,2] < 0.004: # (2) Region wasn't temporary up
                        raise ValueError('1')           # was up, wasn't temporary up and is now low
                        

            elif state_grid[i,j,0] < 0.004:         # (1) Region was low
            
                if state_grid[i,j,1] > 0.004:       # (3) Region is up
                    out_grid[i,j] = 3                   # was low, (became up) and is now up
                    
                elif state_grid[i,j,1] < 0.004:     # (3) Region is now low
                    
                    if state_grid[i,j,2] > 0.004:   # (2) Region was temporary up
                        out_grid[i,j] = 2               # was low, was temporary up and is now low
                    
                    elif state_grid[i,j,2] < 0.004: # (2) Region wasn't temporary up
                        out_grid[i,j] = 1               # was low and stayed permantly low
            else:
                raise ValueError('2')
                            
    return out_grid
    
    
def cleargrid(state_grid):
    [x,y,z]=np.shape(state_grid)
    out_grid=np.zeros((x,y))    
    for i in range(x):
        for j in range(y):
            if state_grid[i,j,1] > 0.004:
                out_grid[i,j] = 3               # permantly up
            elif state_grid[i,j,1] < 0.004:
                if state_grid[i,j,2] > 0.004+state_grid[i,j,0]:
                    out_grid[i,j] = 2.1         # temporary up
                    #print "tempo"
                else:
                    out_grid[i,j] = 1
            else:
                raise ValueError('Error')
                print "ERROR"
                
    return out_grid
# 
def classifygrid(state_grid):
    [x,y,z]=np.shape(state_grid)
    sectiongrid=np.zeros((x,y))
    for i in range(x):
        for j in range(y):
            if state_grid[i,j,2] > 0.004:
                if state_grid[i,j,1] > 0.004:
                    sectiongrid[i,j] = 3.0   #permanent state
                else:
                    sectiongrid[i,j] = 2.1   #temporal state
            elif state_grid[i,j,2] < 0.004:
                sectiongrid[i,j] = 1   #low state
            else:
                raise ValueError('Error')
                print "ERROR"
                
    return sectiongrid
           

###############################################################################
#
# Analysis
#
###############################################################################

# baselevel, 2: endlevel, 3:max()
from matplotlib import cm   
import matplotlib.pyplot as plt

###full
#hirange = np.arange(2,8.1,0.5)*1e-3
#herange = np.arange(7.,0.99,-0.5)*1e-3
#length_range = np.arange(500,1501,10)
#intensity_range = np.arange(250,69,-2)

#full2
hirange = np.arange(10,26,1)*1e-3
herange = np.arange(7.,2.49,-0.25)*1e-3
length_range = np.arange(500,1501,10)
intensity_range = np.arange(250,49,-2)


low_lenge=np.min(length_range)
high_lenge=np.max(length_range)
low_inte=np.min(intensity_range)
high_inte=np.max(intensity_range)




#%%

state_grid=np.zeros((len(herange),len(hirange),len(intensity_range),len(length_range),3))     
i=0
j=0
for he in herange:
    j=0    
    for hi in hirange:
        #state_grid[i,j,:,:,:]=np.load('RUM_Dec_meas_Npp100_HeHi_le500t1500i10msInt50t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000))
        state_grid[i,j,:,:,:]=np.load('RUM_Dec_meas_full2_le500t1500i10msInt50t250i2_He%.2fmV_Hi%.1fmV.npy' %(he*1000,hi*1000))
        print "he:%.2f | hi:%.2f " %(he*1000,hi*1000)
        j=j+1
    i=i+1


glob_low_val_max=np.min(state_grid[:,:,:,:,2])
glob_high_val_max=np.max(state_grid[:,:,:,:,2])

glob_low_val_dete=np.min(state_grid[:,:,:,:,1])
glob_high_val_dete=np.max(state_grid[:,:,:,:,1])


glob_low_val_base=np.min(state_grid[:,:,:,:,0])
glob_high_val_base=np.max(state_grid[:,:,:,:,0])


scale_base=np.zeros((2,2))
scale_base[0]=glob_low_val_base
scale_base[1]=glob_high_val_base
scale_detec=np.zeros((2,2))
scale_detec[0]=glob_low_val_dete
scale_detec[1]=glob_high_val_dete
scale_max=np.zeros((2,2))
scale_max[0]=glob_low_val_max
scale_max[1]=glob_high_val_max
scale_colorcode=np.array([[5,5],[4,4],[3,3],[2,2],[1,1]])


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
plt.figure(1) 
plt.clf()
plt.figure(2) 
plt.clf()          
plt.figure(3) 
plt.clf()
plt.figure(4) 
plt.clf()
plt.figure(5) 
plt.clf()
i=0
j=0
n=1

for he in herange:
    j=0    
    for hi in hirange:
        

        print "he:%.2f | hi:%.2f " %(he*1000,hi*1000)
        
        #        if np.min(state_grid[i,j,:,:,2]) < glob_low_val_max: 
        #            glob_low_val_max=np.min(state_grid[i,j,:,:,2])
        #            print he,hi,glob_low_val_max
        #        if np.max(state_grid[i,j,:,:,2]) > glob_high_val_max: 
        #            glob_high_val_max=np.max(state_grid[i,j,:,:,2])
        #            print he,hi,glob_high_val_max,1
        #    
        #        if np.min(state_grid[i,j,:,:,1]) < glob_low_val_dete: 
        #            glob_low_val_dete=np.min(state_grid[i,j,:,:,1])
        #            print he,hi,glob_low_val_dete
        #        if np.max(state_grid[i,j,:,:,1]) > glob_high_val_dete: 
        #            glob_high_val_dete=np.max(state_grid[i,j,:,:,1])
        #            print he,hi,glob_high_val_dete,1
        #        if np.min(state_grid[i,j,:,:,0]) < glob_low_val_base: 
        #            glob_low_val_base=np.min(state_grid[i,j,:,:,0])
        #            print he,hi,glob_low_val_base
        #        if np.max(state_grid[i,j,:,:,0]) > glob_high_val_base: 
        #            glob_high_val_base=np.max(state_grid[i,j,:,:,0])
        #            print he,hi,glob_high_val_base,1

        #Visualize parameter ranges
        plt.figure(1)
        plt.subplot(len(herange),len(hirange),n)
        plt.text(0.5, 0.5, "'he:%.2f,hi:%.2f,#:%.0f" %(he*1000,hi*1000,n), size=8, rotation=0.,
        ha="center", va="center",
        bbox=dict(boxstyle="square", ec=(1., 1., 1.), fc=(1., 1., 1.),  ))
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])
    
        #Visualize base level
        plt.figure(2)
        plt.subplot(len(herange),len(hirange),n)
        current_grid=state_grid[i,j,:,:,:]  
        current_grid[-1,-1]=glob_low_val_base
        current_grid[-1,-2]=glob_high_val_base
        plt.imshow(np.flipud(current_grid[:,:,0]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none',cmap=cm.Accent)
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])

        #Visualize detector level
        plt.figure(3)
        plt.subplot(len(herange),len(hirange),n)
        current_grid=state_grid[i,j,:,:,:] 
        current_grid[-1,-1]=glob_low_val_dete
        current_grid[-1,-2]=glob_high_val_dete
        plt.imshow(np.flipud(current_grid[:,:,1]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none',cmap=cm.Accent)
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])


        plt.figure(4)
        plt.subplot(len(herange),len(hirange),n)
        current_grid=state_grid[i,j,:,:,:]   
        current_grid[-1,-1]=glob_low_val_max
        current_grid[-1,-2]=glob_high_val_max
        plt.imshow(np.flipud(current_grid[:,:,2]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none',cmap=cm.Accent)
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])
    

        
        # for Visualization of the sectioned values       
        plt.figure(5)
        plt.subplot(len(herange),len(hirange),n)
        current_grid=cleargrid2(state_grid[i,j,:,:,:])   
        current_grid[-1,-1]=5
        current_grid[1,2]=1
        plt.imshow(np.flipud(current_grid[:,:]), aspect='auto', extent = (low_lenge,high_lenge,low_inte,high_inte),interpolation='none',cmap=cm.Accent)
        a=plt.gca()
        a.axes.set_xticklabels([])
        a.axes.set_yticklabels([])
        
        j=j+1
        n=n+1
        if n==3000:
            raise ValueError('dd')

    i=i+1


plt.figure(2)
plt.title('Baselevel') 
plt.tight_layout()
plt.figure(3)
plt.title('Detectorlevel') 
plt.tight_layout()
plt.figure(4)
plt.title('Maximum Values') 
plt.tight_layout()
plt.figure(5)
plt.tight_layout()
#cb=plt.colorbar()
#plt.ylabel('intensity')
#plt.xlabel('length')
print "fertig"        

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

plt.figure(1)
plt.savefig('RUM_Detektor_2popNpp_1_axis.pdf', format='pdf', dpi=1000)
plt.figure(2)
plt.savefig('RUM_Detektor_2popNpp_1_baselevel.pdf', format='pdf', dpi=1000)
plt.figure(3)
plt.savefig('RUM_Detektor_2popNpp_1_detectorlevel.pdf', format='pdf', dpi=1000)
plt.figure(4)
plt.savefig('RUM_Detektor_2popNpp_1_maxlevel.pdf', format='pdf', dpi=1000)
plt.figure(5)
plt.savefig('RUM_Detektor_2popNpp_1_sectioned.pdf', format='pdf', dpi=1000)
plt.figure(14)
plt.savefig('RUM_Detektor_2popNpp_1_baselevel_single.pdf', format='pdf', dpi=1000)
plt.figure(15)
plt.savefig('RUM_Detektor_2popNpp_1_detelevel_single.pdf', format='pdf', dpi=1000)
plt.figure(16)
plt.savefig('RUM_Detektor_2popNpp_1_maxlevel_single.pdf', format='pdf', dpi=1000)
#plt.close('all')
