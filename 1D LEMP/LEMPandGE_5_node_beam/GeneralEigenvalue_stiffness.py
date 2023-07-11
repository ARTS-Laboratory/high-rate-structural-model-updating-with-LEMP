# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 16:17:55 2020

@author: clair
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.linalg as solver
import decimal as dec
import sympy as symp
import json as json
from scipy.signal import butter, lfilter
from scipy.signal import get_window

#%%Set Model Info

######################################################################################
#model info
with open('modelParameters.json') as json_file:
    data = json.load(json_file)
   
for keys, vals in data.items():
    exec(keys + '=vals')

######################################################################################
    
np.set_printoptions(np.set_printoptions(precision=50))
    
# Generating the uniform grid    
beam1 = np.linspace(0, 1, data["nodes"] , endpoint=True)
pin_locations_FEA = beam1
    
    
beam_length = 0.3525     # length of the beam in meters
beam_width = 0.0508   # width of the beam in meters
beam_height = 0.00635 # thickness of the beam in meters
beam_E = 209500000000 # Youngs modules of steel in Pa
beam_density = 7900 # density of steel in kg/m^3
accelerometer_mass = 0.07 # mass of accelerometer in kg
beam_I = (beam_width*beam_height**3)/12 # caclulated moment of inertia
beam_area = beam_width*beam_height
pin_node_rotation_spring = 0 # set the value of the spring at the pinned connection  
    
pin_locations_actual = pin_locations_FEA*beam_length
beam_node_num = pin_locations_FEA.size
beam_element = beam_node_num-1 # calculate the number of elements in the beam
beam_el_length = beam_length/beam_element
 
#%% State 1
matrix_size = (beam_node_num)*2
M = np.zeros((matrix_size,matrix_size))
K = np.zeros((matrix_size,matrix_size))
    
# for each element, add the element matrix into the global matirx
for elem_num in range(0,beam_element):
    if elem_num == (beam_element-1):
        beam_density = beam_density + accelerometer_mass/(beam_area*beam_el_length)
      
       
    # define the mass matrix of a Euler-Bernoulli beam
    M_el = (beam_density*beam_area*beam_el_length)/420* \
    np.matrix([[156,22*beam_el_length,54,-13*beam_el_length], \
                   [22*beam_el_length,4*beam_el_length**2,13*beam_el_length,-3*beam_el_length**2], \
                   [54,13*beam_el_length,156,-22*beam_el_length], \
                   [-13*beam_el_length,-3*beam_el_length**2,-22*beam_el_length,4*beam_el_length**2]])
            
    # define the stiffness matrix of a Euler-Bernoulli beam
    K_el = (beam_E*beam_I)/beam_el_length**3* \
    np.matrix([[12,6*beam_el_length,-12,6*beam_el_length], \
                   [6*beam_el_length,4*beam_el_length**2,-6*beam_el_length,2*beam_el_length**2], \
                   [-12,-6*beam_el_length,12,-6*beam_el_length], \
                   [6*beam_el_length,2*beam_el_length**2,-6*beam_el_length,4*beam_el_length**2]])
                
    n = (elem_num)*2
    M[n:n+4,n:n+4] = np.add(M[n:n+4,n:n+4],M_el)
    K[n:n+4,n:n+4] = np.add(K[n:n+4,n:n+4],K_el)    
            
# for the fixed end on the left side, u_1 and u_2 = 0, so we can remove these columns
# and rows form the matrixes. 
# apply the boundary conditions
K[0,0]=1e10
K[1,1]=1e10
# K[2,2]=data["pin_equivalent_stiffness"]
K[7,7]=data["pin_equivalent_stiffness"]
# K[6,6]=data["pin_equivalent_stiffness"]
# K[8,8]=data["pin_equivalent_stiffness"]
# K[10,10]=data["pin_equivalent_stiffness"]


#%% State 2 Function
   
def State2(pin_loc, number_modes, number_nodes):
    spring_loc= np.int(np.round((pin_loc*(number_nodes-1))))
    
#################################### STATE 2 ########################################
    del_K = np.zeros((matrix_size,matrix_size))
    del_K.itemset(((spring_loc)*2, (spring_loc)*2), data["pin_stiffness_change"])
    
    mM2= M
    kK2= K + del_K

    # Calculation of the natural frequencies. 
    eigvals,eigvects = sp.linalg.eigh(kK2,mM2)
    eigvals=np.expand_dims(np.real(eigvals), axis=0)
    FEA_wn = np.sort(np.real(np.squeeze(np.sqrt(eigvals)))) # Natural frequencies, rad/s
    Frequencies = FEA_wn/(2*np.pi) # Natural freq in Hz

    
    return(Frequencies[0:number_modes])
