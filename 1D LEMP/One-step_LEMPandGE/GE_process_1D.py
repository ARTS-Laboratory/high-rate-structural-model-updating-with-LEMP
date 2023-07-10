 # -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 14:44:15 2021

@author: OGUNNIYI
"""
'''
This is code uses the generalized eigenvalue solution to solve for a single change on a beam
'''

import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy.solvers import solve
from sympy import Symbol
import scipy.linalg as solve
import sympy as sp
import math
import time
import random

import matplotlib.pyplot as plt
import scipy as sp
import scipy.linalg as solver
import decimal as dec
import sympy as symp
import json as json
from scipy.signal import butter, lfilter
from scipy.signal import get_window


from datetime import datetime
import secular_equation_functionV4 as solve_sec

''' code written to work with 4 element beam (5 nodes) for first five modes
'''

#%% define beam properties

node = 5 # beam number of nodes
beam_element_num = node-1 # number of elements 
DOF = node*2 # degree of freedom for the system
mode = 5


#%% Material properties of beam

Density = 7900
d = Density
Youngs_Modulus = 209500000000
Y = Youngs_Modulus
beam_length = 0.3525     # length of the beam in meters
beam_width = 0.0508   # width of the beam in meters
beam_height = 0.00635 # thickness of the beam in meters
accelerometer_mass = 0.07 # mass of accelerometer in kg
beam_I = (beam_width*beam_height**3)/12 # caclulated moment of inertia
beam_area = beam_width*beam_height
pin_node_rotation_spring = 0 # set the value of the spring at the pinned connection

Cross_sectional_Area = beam_area
A = Cross_sectional_Area
Total_length = beam_length
TL = Total_length
Element_length = beam_length/beam_element_num
EL = Element_length
pstif = pin_equivalent_stiffness = 10000000000

#%% Building up the global mass and stiffness matrix

matrix_size = (node)*2
M1 = np.zeros((matrix_size,matrix_size))
K1 = np.zeros((matrix_size,matrix_size))

for i in range(0,beam_element_num):
    if i == (beam_element_num-1):
        d = d + accelerometer_mass/(beam_area*EL)
        
# define the mass matrix of a Euler-Bernoulli beam
    M_el = (d*A*EL)/420* \
    np.matrix([[156,22*EL,54,-13*EL], \
                   [22*EL,4*EL**2,13*EL,-3*EL**2], \
                   [54,13*EL,156,-22*EL], \
                   [-13*EL,-3*EL**2,-22*EL,4*EL**2]])
            
    # define the stiffness matrix of a Euler-Bernoulli beam
    K_el = (Y*beam_I)/EL**3* \
    np.matrix([[12,6*EL,-12,6*EL], \
                   [6*EL,4*EL**2,-6*EL,2*EL**2], \
                   [-12,-6*EL,12,-6*EL], \
                   [6*EL,2*EL**2,-6*EL,4*EL**2]])
    
    n = (i)*2
    M1[n:n+4,n:n+4] = np.add(M1[n:n+4,n:n+4],M_el)
    K1[n:n+4,n:n+4] = np.add(K1[n:n+4,n:n+4],K_el)  
       
# for the fixed end on the left side, u_1 and u_2 = 0, so we can remove these columns
# and rows form the matrixes. 
# apply the boundary conditions

K1[0,0]=1e10
K1[1,1]=1e10 
    

#%% Initial state values for eigenvalues and eigenvectors

def eig_soln(M1,K1):
    # fuction to calculation of the natural frequencies using the generalized eigenvalue approach
    #also returns matrix of corresponding eigenvectors U1
    eigvals,eigvects = solve.eigh(K1,M1)
    eigvals=np.expand_dims(np.real(eigvals), axis=0)
    
    U1=eigvects
    wn_squared= eigvals
    
    FEA_wn= np.sqrt(wn_squared)
    Frequencies = FEA_wn/(2*np.pi) # Natural freq in Hz
    
    return(wn_squared,U1, Frequencies)

state1_values= eig_soln(M1,K1)
wn_squared= state1_values[0]
old_frequencies= state1_values[2]
U1= state1_values[1]

old_frequencies = old_frequencies[:,0:mode]
# print(old_frequencies)


#%% GE process
# ALTERED STATE

# Adding stiffness at node 4, DOF 7 (deflection), then del_K will be

spring_node = 3 
del_K = np.zeros((matrix_size,matrix_size))
del_K.itemset((((spring_node)*2) , ((spring_node)*2)), 1e10) # change in stiffness
# print(del_K)

mM2= M1 # new mass matrix
kK2= K1 + del_K # new stiffness matrix

# Calculation of the natural frequencies. 
eigvals,eigvects = sp.linalg.eigh(kK2,mM2)
eigvals=np.expand_dims(np.real(eigvals), axis=0)
FEA_wn = np.sort(np.real(np.squeeze(np.sqrt(eigvals)))) # Natural frequencies, rad/s
Frequencies = FEA_wn/(2*np.pi) # Natural freq in Hz

print(Frequencies[0:5]) # print first five mode frequencies

