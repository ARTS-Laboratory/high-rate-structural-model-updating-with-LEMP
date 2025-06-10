# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:54:08 2020

@author: clair
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import json as json
from datetime import datetime
import math

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 17})
plt.rcParams.update({'axes.titlesize': 17})
plt.rcParams.update({'axes.labelsize': 17})
plt.rcParams.update({'lines.markersize': 17})
plt.rcParams.update({'legend.fontsize': 17})
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 



'''
Code to compare state estimation using the general eigenvalue and LEMP algorithm
'''

#%%
#MODEL INFO
nodes= 5
modes= 6

pin_equivalent_stiffness=1e10
pin_stiffness_change=-7e9

modelParameters = {'nodes':nodes,
        'modes':modes, 
        'pin_equivalent_stiffness': pin_equivalent_stiffness,
        'pin_stiffness_change': pin_stiffness_change
        }

with open('modelParameters.json', 'w') as outfile:
    json.dump(modelParameters, outfile)

#%%Import codes
import GeneralEigenvalue_stiffness as GES
import LEMP_simplified_procedure as lemp2

#%%
#State 1 Info
wn_squared= lemp2.wn_squared
#wn_squared=wn_squared[0,0:modes]
old_freq= lemp2.old_frequencies
old_freq=old_freq[0,0:modes]

U1=lemp2.U1

# print(wn_squared)
print(old_freq)
#%%
#Test LEMP vs general eigen solution at midpoint of beam
test_freq_GES=GES.State2(pin_loc=0.5, number_modes=modes, number_nodes=nodes)
test_freq_LEMP=lemp2.LEMP(pin_loc=0.5, number_modes=modes, number_nodes=nodes)



#TEST using 0-1 scale for beam locations
#Test modal shapes acording to LEMP as roller position changes
#plots 1-5 LEMP, plots 6-10 Eigensolver
pin_locations_FEA = np.linspace(0, 1, nodes , endpoint=True)


pin_locations= pin_locations_FEA[1:]


new_frequencies = np.zeros((modes,nodes-1))
new_frequencies_GE = np.zeros((modes,nodes-1))


for i in (pin_locations):
    new_freq= lemp2.LEMP(pin_loc=i, number_modes=modes, number_nodes=nodes)
    new_freq_GE= GES.State2(pin_loc=i, number_modes=modes, number_nodes=nodes)
    
    index = np.where(pin_locations==i)
    #index = np.int(indexes[0])
    
    for j in range (modes-1):
        new_frequencies[j, index]= new_freq[j]
        new_frequencies_GE[j, index]= new_freq_GE[j]
     

print(new_frequencies)
print(new_frequencies_GE)
#%% calculating MEA

error = new_frequencies - new_frequencies_GE

print(error)
summation_0= sum(error[0,:])
summation_1= sum(error[1,:])
summation_2= sum(error[2,:])
summation_3= sum(error[3,:])
summation_4= sum(error[4,:])


MEA_0 = summation_0/len(error[0,:])
MEA_1 = summation_1/len(error[1,:])
MEA_2 = summation_2/len(error[2,:])
MEA_3 = summation_3/len(error[3,:])
MEA_4 = summation_4/len(error[4,:])

  
#%% calculating signal to noise ratio
SNR = np.zeros(5)

SNR[0] = math.log10(abs(np.mean(new_frequencies_GE[0,:])/np.mean(error[0,:])))*10
SNR[1] = math.log10(abs(np.mean(new_frequencies_GE[1,:])/np.mean(error[1,:])))*10
SNR[2] = math.log10(abs(np.mean(new_frequencies_GE[2,:])/np.mean(error[2,:])))*10
SNR[3] = math.log10(abs(np.mean(new_frequencies_GE[3,:])/np.mean(error[3,:])))*10
SNR[4] = math.log10(abs(np.mean(new_frequencies_GE[4,:])/np.mean(error[4,:])))*10

r=abs(np.mean(new_frequencies_GE[0,:])/np.mean(error[0,:]))
print(r)

mode_ = np.array([1,2,3,4,5])

plt.figure(figsize=(6,4)) 
plt.plot(mode_,SNR,'-o',markersize= 15,label='mode 1')   
plt.locator_params(axis="x", nbins=5)
plt.ylim(12,35)
plt.xlabel('mode number')
plt.ylabel('SNR$_{dB}$')
plt.legend(framealpha=1,fontsize=17,loc='upper right', fancybox=True,ncol=1, bbox_to_anchor=(1,1))
plt.tight_layout()
plt.savefig('SNR1.pdf', dpi=500)
plt.show()


# #%% Error plot for all mode

zero_error = np.zeros(len(error[0,:]))

per= (error/new_frequencies_GE)*100


plt.figure(figsize=(6,4))    
plt.plot(pin_locations_FEA[1:], per[0,:], '-d',markersize= 8,label='mode 1')
plt.plot(pin_locations_FEA[1:], per[1,:], '-*',markersize= 8,label='mode 2')
plt.plot(pin_locations_FEA[1:], per[2,:], '-v',markersize= 8,label='mode 3')
plt.plot(pin_locations_FEA[1:], per[3,:], '-s',markersize= 8,label='mode 4')
plt.plot(pin_locations_FEA[1:], per[4,:], '-o',markersize= 8,label='mode 5')
plt.locator_params(axis="x", nbins=6)
plt.locator_params(axis="y", nbins=8)
plt.xlabel('pin positions % of beam')
plt.ylabel('frequency error (%)')
plt.legend(framealpha=1,fontsize=17,loc='upper right', fancybox=True,ncol=2, bbox_to_anchor=(1,1))
plt.tight_layout()
plt.savefig('error plot.pdf', dpi=500)
plt.show()


#%% Plot the results         
  
plt.figure(figsize=(6.9,2))    
plt.plot(pin_locations_FEA[0], old_freq[0], 'bo',markersize=10, label='fixed end')
plt.plot(pin_locations_FEA[1:], new_frequencies_GE[0,:], '-o',markersize=10,label='Generalized Eigenvalue')
plt.plot(pin_locations_FEA[1:], new_frequencies[0,:], '-s',markersize=6, label='LEMP')
# plt.xlabel('pin positions % of beam')
plt.ylabel('frequency (Hz)')
# plt.legend(framealpha=1, fontsize=20, loc=2)
plt.tight_layout()
plt.savefig('mode 1.pdf', dpi=500)
plt.show()


plt.figure(figsize=(6.7,2))
plt.plot(pin_locations_FEA[0], old_freq[1], 'bo',markersize=10, label='original state')
plt.plot(pin_locations_FEA[1:], new_frequencies_GE[1,:], '-o',markersize=10,label='GE')
plt.plot(pin_locations_FEA[1:], new_frequencies[1,:], '-s',markersize=6,label='LEMP')
# plt.xlabel('pin positions % of beam')
# plt.ylabel('frequency (Hz)')
# plt.legend(framealpha=1, fontsize=20)
plt.tight_layout()
plt.savefig('mode 2.pdf', dpi=500)
plt.show()


plt.figure(figsize=(7,2))
plt.plot(pin_locations_FEA[0], old_freq[2], 'bo',markersize=10, label='original state')
plt.plot(pin_locations_FEA[1:], new_frequencies_GE[2,:], '-o',markersize=10,label='GE')
plt.plot(pin_locations_FEA[1:], new_frequencies[2,:], '-s',markersize=6,label='LEMP')

# plt.xlabel('pin positions % of beam')
plt.ylabel('frequency (Hz)')
# plt.legend(framealpha=1, fontsize=20)
plt.tight_layout()
plt.savefig('mode 3.pdf', dpi=500)
plt.show()


plt.figure(figsize=(6.8,2.3))
plt.plot(pin_locations_FEA[0], old_freq[3], 'bo',markersize=10, label='original state')
plt.plot(pin_locations_FEA[1:], new_frequencies_GE[3,:], '-o',markersize=10,label='GE')
plt.plot(pin_locations_FEA[1:], new_frequencies[3,:], '-s',markersize=6,label='LEMP')
plt.xlabel('pin positions % of beam')
# plt.ylabel('frequency (Hz)')
# plt.legend(framealpha=1, fontsize=20)
plt.tight_layout()
plt.savefig('mode 4.pdf', dpi=500)
plt.show()


plt.figure(figsize=(7,2.3))
plt.plot(pin_locations_FEA[0], old_freq[4], 'bo',markersize=10, label='original state')
plt.plot(pin_locations_FEA[1:], new_frequencies_GE[4,:], '-o',markersize=10,label='GE')
plt.plot(pin_locations_FEA[1:], new_frequencies[4,:], '-s',markersize=6,label='LEMP')
plt.xlabel('pin positions % of beam')
plt.ylabel('frequency (Hz)')
# plt.legend(framealpha=1, fontsize=20)
plt.tight_layout()
plt.savefig('mode 5.pdf', dpi=500)
plt.show()
 
 
