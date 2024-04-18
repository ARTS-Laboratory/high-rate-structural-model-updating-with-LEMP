# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:16:36 2022

@author: PXIe-8880
"""

import numpy as np
import scipy.linalg as linalg


Stiff = np.loadtxt('stiff.lvm',delimiter=',')
stfl = int(len(Stiff)**.5)
Stiff = np.reshape(Stiff,(stfl,stfl))
Mass = np.loadtxt('mass.lvm',delimiter=',') 
Mass = np.reshape(Mass,(stfl,stfl))

eigenvals, eigenvectors = linalg.eigh(Stiff, Mass)


eigvalout=''
for val in eigenvals:
    eigvalout += str(val)+','
eigvalout=eigvalout[:-1]
with open('eigenvals.txt','w') as file:
    file.write(eigvalout)
    

eigvectout=''
for vect in eigenvectors.ravel():
    eigvectout += str(vect)+','
eigvectout=eigvectout[:-1]
with open('eigenvects.txt','w') as file:
    file.write(eigvectout)