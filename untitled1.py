#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 14:34:10 2019

@author: shilu
"""

    if u[0].shape[0] == 3:
        
        for sweeps in range(s2):
            u = Rich(u, rhs, alpha)
             # u =Jacobi(u, rhs, alpha, L_list[count]) 
            
    else:
        
        for sweeps in range(s1):
            
            #smoother_plot(u)        

            u = Rich(u, rhs, alpha)
            
            #smoother_plot(u)
             #u = Jacobi(u, rhs, alpha, L_list[count])
             
  
    if u[0].shape[0] != 3:
        
        #print 'dasd'

        
        rhs = residue_sqrt_alpha(rhs, u, alpha)
        
        rhs = Res_Injection(rhs)
        

        uc = np.zeros((4, rhs[0].shape[0], rhs[0].shape[0]))
        

        
        uc = VCycle(uc, rhs, s1, s2, alpha,  count + 1)
        
        u = u + Interpolation(uc)
            

    return u
