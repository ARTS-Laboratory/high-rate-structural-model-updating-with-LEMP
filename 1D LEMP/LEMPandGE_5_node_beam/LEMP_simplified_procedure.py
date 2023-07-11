# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 09:08:57 2020

@author: clair
"""

import numpy as np
import matplotlib.pyplot as plt
import decimal as dec
import sympy as sp
import scipy.linalg as solve
import scipy as sci
import math
import json as json

#%%Set info from json files
np.set_printoptions(np.set_printoptions(precision=50))  

#model info
with open('modelParameters.json') as json_file:
    data = json.load(json_file)
   
for keys, vals in data.items():
    exec(keys + '=vals')

#%%Functions
def state2(spring_node, modes_number, nodes_number):
    
    del_K=np.zeros((nodes_number*2,nodes_number*2))
    del_K.itemset(((spring_node*2), (spring_node*2)), data["pin_stiffness_change"])

    return(del_K)
    
def eig_soln(M,K):
    # Calculation of the natural frequencies using the generalized eigenvalue approach
    #also returns matrix of corresponding eigenvectors U1
    eigvals,eigvects = solve.eigh(K,M)
    eigvals=np.expand_dims(np.real(eigvals), axis=0)
    
    U1=eigvects
    wn_squared= eigvals
    
    FEA_wn= np.sqrt(wn_squared)
    Frequencies = FEA_wn/(2*np.pi) # Natural freq in Hz
    
    return(wn_squared,U1, Frequencies)
 

#%%State 1 info
######################################################################################
    
#np.set_printoptions(np.set_printoptions(precision=50))
    
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

#UNDEFORMED STATE- will always be the same 
#so no need for function here
#save as json file for easy reference

state1_measurements= eig_soln(M,K)
wn_squared= state1_measurements[0]
old_frequencies= state1_measurements[2]
U1= state1_measurements[1]
    
###############################################################################
    
#Make this a function!
#DEFORMED STATE
def LEMP(pin_loc, number_modes, number_nodes):
    
    spring_loc= np.int(np.round((pin_loc*(number_nodes-1))))
    
    #DEFORMED STATE    
    del_K= state2(spring_node=spring_loc, modes_number=number_modes, nodes_number=number_nodes)
    
    #Spectral Decomposition of del_K
    #alpha values and tie vectors (T)
    alpha, T = np.linalg.eig(del_K)
    
    #only connecting nodes contribute to state changes so modifications necessary...
    v_matrix= U1.T @ T
    v_vect_contribute= v_matrix[:,((spring_loc*2))]

    #define how many modes to include i.e. set truncation
    v_vect=v_vect_contribute[0:number_modes]
    wn_squared_state1=wn_squared[:,0:number_modes]
    
    
    #creates an equation using summation
    #sp.init_printing(use_unicode=True)
    
    summation = 0
    vr, wr, o2 = sp.symbols('vr wr o2')
    expr = ((vr * vr)/ (wr-o2))
            
    for r in range (number_modes):
        v_r= v_vect[r]
        wn_r= wn_squared_state1[0,r]
                
        part = expr.subs(vr, v_r).subs(wr, wn_r)
        print(part)
        
        summation= summation + part
        print(summation)
    
        eq= sp.Eq(summation, -1/data["pin_equivalent_stiffness"])
        print(eq)
           
    sol= sp.solveset(eq, o2)
    print(sol)  
        
    modes= list(sol)
    
    new_omegas_squared= np.zeros(len(modes))
    new_omegas= np.zeros(len(modes))
    new_frequencies= np.zeros(len(modes))
    
    for k in range (len(modes)):
        new_omegas_squared[k]= modes[k]
        new_omegas[k]= (np.sqrt(new_omegas_squared[k]))
        new_frequencies[k] = new_omegas[k]/(2*np.pi) # Natural freq in Hz
    new_frequencies.sort()

    return(new_frequencies)
    
