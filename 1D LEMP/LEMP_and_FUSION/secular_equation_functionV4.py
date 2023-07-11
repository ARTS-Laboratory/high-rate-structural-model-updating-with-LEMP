
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:54:35 2021

@author: eaogu
"""

import numpy as np

def secular_eq_solver(mode,d,v,p):
    sol=np.zeros(mode)
    for i in range(1,mode+1):
        k = i
        dk = d[0,k-1]
        
        
        if k<mode:
            dk1 = d[0,k]
            lamda = (dk+dk1)/2
            # print(lamda)
            'To check if f_y is positive or negative'
            f_ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-lamda)
                # print(q)
            f_y = f_ + p
            
            'Now, applying equation 40 and 41, choosing our initial guess'
          
            g_x= f_y-v[k-1]/(d[0,k-1]-lamda)
            
            Tr = d[0,k] - d[0,k-1]
            c1 = g_x
            
            'first a,b and c to obtain y initial guess' 'equation 43 and 44'
            '''
            HERE
            IS 
            THE 
            CHANGE
            '''
            # a1 = c1*Tr + (v[k-1]+v[k])
            # b1 = v[k-1]*Tr 
            
            if f_y >= 0:
                a1 = c1*Tr + (v[k-1]+v[k])
                b1 = v[k-1]*Tr
            elif f_y < 0:
                a1 = -c1*Tr+ (v[k-1]+v[k])
                b1 = -v[k]*Tr
            
            if a1 <= 0:
                Tao = (a1-np.sqrt(a1**2-4*b1*c1))/2*c1
            elif a1 > 0:
                Tao = 2*b1/(a1 + np.sqrt(a1**2-4*b1*c1))
            '''here is the changes'''
            '''
            .
            .
            .
            .
            '''
            # init_y = Tao + d[0,k-1]
            # init_y = Tao + d[0,k]
            if f_y>=0:
                init_y = Tao + d[0,k-1]
                
            elif f_y<0:
                init_y = Tao + d[0,k]
            
            'computing correction n to y for better approximation'
            Tr_k = d[0,k-1] - init_y
            Tr_k_1 = d[0,k] - init_y
            
            f_ = 0
            f__ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-init_y)
                f__ = f__ + v[q]/((d[0,q]-init_y)**2)
            f_y = f_ + p
            f_y_ = f__
            
            o_ = 0
            for q in range(k,mode):
                o_ = o_ + v[q]/((d[0,q]-init_y)**2)
                # print(q)
            o_y_= o_
            u_ = 0
            for q in range(mode):
                # print(q)
                u_ = u_ + v[q]/((d[0,q]-init_y)**2)
            u_y_= u_
            'second a, b and c for final correction'
            
            a2 = (Tr_k + Tr_k_1)*f_y - Tr_k*Tr_k_1*f_y_
            b2 = Tr_k*Tr_k_1*f_y
            
            # c2 = f_y - Tr_k*u_y_ - Tr_k_1*o_y_
            if f_y>=0:
                c2 = f_y - Tr_k_1*f_y_ - u_y_*(d[0,k-1]-d[0,k])
            elif f_y<0:
                c2 = f_y - Tr_k*f_y_ - o_y_*(d[0,k]-d[0,k-1])
            
            
            if a2 <= 0:
                n_cor = (a2-np.sqrt(a2**2-4*b2*c2))/2*c2
            elif a2 > 0:
                n_cor = 2*b2/(a2 + np.sqrt(a2**2-4*b2*c2))
                
            lamda_new = init_y + n_cor
            
            
            'computing second correction n to y for better approximation'
            init_y = lamda_new
            
            Tr_k = d[0,k-1] - init_y
            Tr_k_1 = d[0,k] - init_y
            
            f_ = 0
            f__ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-init_y)
                f__ = f__ + v[q]/((d[0,q]-init_y)**2)
            f_y = f_ + p
            f_y_ = f__
            
            o_ = 0
            for q in range(k,mode):
                o_ = o_ + v[q]/((d[0,q]-init_y)**2)
                # print(q)
            o_y_= o_
            u_ = 0
            for q in range(mode):
                # print(q)
                u_ = u_ + v[q]/((d[0,q]-init_y)**2)
            u_y_= u_
            'second a, b and c for final correction'
            
            a2 = (Tr_k + Tr_k_1)*f_y - Tr_k*Tr_k_1*f_y_
            b2 = Tr_k*Tr_k_1*f_y
            
            # c2 = f_y - Tr_k*u_y_ - Tr_k_1*o_y_
            if f_y>=0:
                c2 = f_y - Tr_k_1*f_y_ - u_y_*(d[0,k-1]-d[0,k])
            elif f_y<0:
                c2 = f_y - Tr_k*f_y_ - o_y_*(d[0,k]-d[0,k-1])
            
            
            if a2 <= 0:
                n_cor = (a2-np.sqrt(a2**2-4*b2*c2))/2*c2
            elif a2 > 0:
                n_cor = 2*b2/(a2 + np.sqrt(a2**2-4*b2*c2))
                
            lamda_new = init_y + n_cor
            
            'computing third correction n to y for better approximation'
            init_y = lamda_new
            
            Tr_k = d[0,k-1] - init_y
            Tr_k_1 = d[0,k] - init_y
            
            f_ = 0
            f__ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-init_y)
                f__ = f__ + v[q]/((d[0,q]-init_y)**2)
            f_y = f_ + p
            f_y_ = f__
            
            o_ = 0
            for q in range(k,mode):
                o_ = o_ + v[q]/((d[0,q]-init_y)**2)
                # print(q)
            o_y_= o_
            u_ = 0
            for q in range(mode):
                # print(q)
                u_ = u_ + v[q]/((d[0,q]-init_y)**2)
            u_y_= u_
            'second a, b and c for final correction'
            
            a2 = (Tr_k + Tr_k_1)*f_y - Tr_k*Tr_k_1*f_y_
            b2 = Tr_k*Tr_k_1*f_y
            
            # c2 = f_y - Tr_k*u_y_ - Tr_k_1*o_y_
            if f_y>=0:
                c2 = f_y - Tr_k_1*f_y_ - u_y_*(d[0,k-1]-d[0,k])
            elif f_y<0:
                c2 = f_y - Tr_k*f_y_ - o_y_*(d[0,k]-d[0,k-1])
            
            
            if a2 <= 0:
                n_cor = (a2-np.sqrt(a2**2-4*b2*c2))/2*c2
            elif a2 > 0:
                n_cor = 2*b2/(a2 + np.sqrt(a2**2-4*b2*c2))
                
            lamda_new = init_y + n_cor
            
            sol[i-1] = lamda_new
        
        elif k==mode:
            t = 0
            for i in range(mode):
                t = t + v[i]**2
                
            d6 = d[0,k-1] + t/p
            
            lamda = (d6 + d[0,k-1])/2
            # print(lamda)
            
            'To check if f_y is positive or negative'
            f_ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-lamda)
            f_y = f_ + p
            # print(f_y)
            'Now, applying equation 40 and 41, choosing our initial guess'
            # g_ = 0
            # for q in range(0,mode,1):
            #     if q != k-1:
            #         g_ = g_ + v[q]/(d[0,q]-lamda)
            #         # print(q)        
            # g_x = g_ + p
            
            g_x= f_y-v[k-1]/(d[0,k-1]-lamda)
        
            
            Tr = d[0,k-1] - d[0,k-2]
            c1 = g_x
            
            'first a,b and c to obtain y initail guess' 'equation 43 and 44'
            if f_y <= 0:
                a1 = c1*Tr + (v[k-2]+v[k-1])
                b1 = v[k-1]*Tr
            elif f_y > 0:
                a1 = -c1*Tr+ (v[k-2]+v[k-1])
                b1 = -v[k-1]*Tr
            
            # if a1 <= 0:
            #     Tao = (a1-np.sqrt(a1**2-4*b1*c1))
            # elif a1 > 0:
            #     Tao = 2*b1/(a1 + np.sqrt(a1**2-4*b1*c1))
            
            # a1 = -c1*Tr+ (v[k-2]+v[k-1])
            # b1 = -v[k-1]*Tr
            # print(a1,b1,c1)
            if a1 <= 0:
                Tao = (a1+np.sqrt(a1**2-4*b1*c1))/2*c1
            elif a1 > 0:
                Tao = 2*b1/(a1 - np.sqrt(a1**2-4*b1*c1))
          
            # print(Tao)
            init_y = Tao + d[0,k-1]
            # print(init_y)
            'computing correction n to y for better approximation'
            Tr_k = d[0,k-1] - init_y
            Tr_k_1 = d6 - init_y
            
            f_ = 0
            f__ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-init_y)
                f__ = f__ + v[q]/((d[0,q]-init_y)**2)
            f_y = f_ + p
            f_y_ = f__
            
            u_ = 0
            for q in range(mode):
                # print(q)
                u_ = u_ + v[q]/((d[0,q]-init_y)**2)
            u_y_= u_
            'second a, b and c for final correction'
            
            a2 = (Tr_k + Tr_k_1)*f_y - Tr_k*Tr_k_1*f_y_
            b2 = Tr_k*Tr_k_1*f_y
            
            c2 = f_y - Tr_k*u_y_ - (v[k-1]/Tr_k_1)
            
            
            if a2 <= 0:
                n_cor = (a2-np.sqrt(a2**2-4*b2*c2))/2*c2
            elif a2 > 0:
                n_cor = 2*b2/(a2 + np.sqrt(a2**2-4*b2*c2))
                
            lamda_new = init_y + n_cor
            
            'computing second correction n to y for better approximation'
            init_y = lamda_new
            
            Tr_k = d[0,k-1] - init_y
            Tr_k_1 = d6 - init_y
            
            f_ = 0
            f__ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-init_y)
                f__ = f__ + v[q]/((d[0,q]-init_y)**2)
            f_y = f_ + p
            f_y_ = f__
            
            u_ = 0
            for q in range(mode):
                # print(q)
                u_ = u_ + v[q]/((d[0,q]-init_y)**2)
            u_y_= u_
            'second a, b and c for final correction'
            
            a2 = (Tr_k + Tr_k_1)*f_y - Tr_k*Tr_k_1*f_y_
            b2 = Tr_k*Tr_k_1*f_y
            
            c2 = f_y - Tr_k*u_y_ - (v[k-1]/Tr_k_1)
            
            
            if a2 <= 0:
                n_cor = (a2-np.sqrt(a2**2-4*b2*c2))/2*c2
            elif a2 > 0:
                n_cor = 2*b2/(a2 + np.sqrt(a2**2-4*b2*c2))
                
            lamda_new = init_y + n_cor
            
            'computing third correction n to y for better approximation'
            init_y = lamda_new
            
            Tr_k = d[0,k-1] - init_y
            Tr_k_1 = d6 - init_y
            
            f_ = 0
            f__ = 0
            for q in range(mode):
                f_ = f_ + v[q]/(d[0,q]-init_y)
                f__ = f__ + v[q]/((d[0,q]-init_y)**2)
            f_y = f_ + p
            f_y_ = f__
            
            u_ = 0
            for q in range(mode):
                # print(q)
                u_ = u_ + v[q]/((d[0,q]-init_y)**2)
            u_y_= u_
            'second a, b and c for final correction'
            
            a2 = (Tr_k + Tr_k_1)*f_y - Tr_k*Tr_k_1*f_y_
            b2 = Tr_k*Tr_k_1*f_y
            
            c2 = f_y - Tr_k*u_y_ - (v[k-1]/Tr_k_1)
            
            
            if a2 <= 0:
                n_cor = (a2-np.sqrt(a2**2-4*b2*c2))/2*c2
            elif a2 > 0:
                n_cor = 2*b2/(a2 + np.sqrt(a2**2-4*b2*c2))
                
            lamda_new = init_y + n_cor
            
            sol[mode-1]=lamda_new
    
    return sol
