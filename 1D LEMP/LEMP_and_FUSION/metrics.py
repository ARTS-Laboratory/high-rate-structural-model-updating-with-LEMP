# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:04:32 2023

@author: AVEREEN
"""

import numpy as np
import scipy as sp

def clean_nans(s):
    for i,y in enumerate(s):
        if np.isnan(y):
            s[i]=s[i-1]
    return s

def resample(y,a,b):
    return clean_nans(sp.signal.resample_poly(y,a,b))

def one_to_one(y1,y2):
    up, down = min([len(y1),len(y2)]), max([len(y1),len(y2)])
    inputs = {str(len(y1)):y1,str(len(y2)):y2}
    ys1 = clean_nans(sp.signal.resample_poly(inputs[str(down)],up,down))
    ys2 = inputs[str(up)]
    return ys1, ys2

def TRAC(y1,y2):
    ys1, ys2 = one_to_one(y1,y2)
    s1, s2 = np.matrix(ys1), np.matrix(ys2)
    return np.sum(((s1@s2.T)**2)/((s1@s1.T)*(s2@s2.T)))

def MSE(y1,y2):
    ys1, ys2 = one_to_one(y1,y2)
    s1, s2 = np.matrix(ys1), np.matrix(ys2)
    return np.mean(np.square(s1-s2))

def SNR(y1,y2):
    ys1, ys2 = one_to_one(y1,y2)
    s1, s2 = np.matrix(ys1), np.matrix(ys2)
    return 10*np.log10(np.sum(np.square(s2))/np.sum(np.square(np.subtract(s1,s2))))