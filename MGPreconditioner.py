#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:18:10 2019

@author: shilu
"""


from Transfer import Interpolation, Res_Injection, FW_Restriction
from Residue import residue_sqrt_alpha
from PlotRountines import smoother_plot, convergence_plot
import numpy as np
from functions import Exy, Xexy, Yexy, XYexy, sin_soln,Linear2x, Xlinear2x, Ylinear2x, Zero
from PreProcess import pre_processing
from Triangle_Matrices import h_boundary, A_connection, Amatrix, Lmatrix
    

def VCycle(u,rhs,  s1, s2, alpha, num_data):
    """ This routine implements the recurssion version of V-cycle
        input: current approximation of u
               rhs 
               s1: number of iterations for relaxtions applied downwards
               s2: number of iterations for relaxtions applied upwards
               alpha
        outpu: new approximation of u
    """
    
    from MGSmoother import Rich
    from UzawaSmoother import Uzawa
    

    
    
    
    
    if u[0].shape[0] != 3:
        
        l = int(np.log2(u[0].shape[0]-1))
        
        print l
        grid = pre_processing(l,num_data)
        Ama = Amatrix(grid, l, num_data)[0]
        Lma = Lmatrix(grid, l, num_data)
        del grid
        
        for sweeps in range(s1):
            
            smoother_plot(u[0])
#            print u,'before'
            u = Uzawa(u, rhs, alpha, Ama, Lma)
#            u = Rich(u, rhs,alpha, Ama)
#            print u,'after'
#            smoother_plot(u[0])   
        
        rhs1 = residue_sqrt_alpha(rhs, u, alpha, Ama)

        rhs1 = FW_Restriction(rhs1)
        
        uc = np.zeros((4, rhs1[0].shape[0], rhs1[0].shape[0]))
    
        h = 1/float(rhs1[0].shape[0]-1)
        
        x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
        

        uc = VCycle(uc, rhs1, s1, s2, alpha, num_data)
        
        print 'up'
        u = u + Interpolation(uc)
        
        #POST SMOOTHING???
        
    else:
        
    
        grid = pre_processing(1, num_data)
        Ama = Amatrix(grid, 1, num_data)[0]
        Lma = Lmatrix(grid, 2, num_data)
        for sweeps in range(s2):
             
             u = Uzawa(u, rhs, alpha, Ama, Lma)
#            u = Rich(u, rhs,alpha, Ama)
#             print u, 'bottom'
#             print u[0]
#             smoother_plot(u[])
             
        
    return u



def MG(i, alpha, num_data, num_cyc):

    
    
    n = 2**i +1
    int_n = n-2
    
    h=1/float(n-1)
    
    c = Exy
    
    g1 = Xexy
    
    g2 = Yexy
    
    w = XYexy
    
   
    x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    
    
    grid = pre_processing(i, num_data)
    
    A_global = A_connection(grid,i, num_data)
    
    print 'done A_global'
    
    Ama, d= Amatrix(grid, i, num_data)
    
    print 'done d'

    b = h_boundary(grid, A_global, d, i,num_data, c, g1, g2, w,alpha)
    
#    bb = h_boundary(grid, A_global, d, i,num_data, c, g1, g2, w,1e-8)
    
#    b = h_boundary(grid, A_global, d, i,num_data, Exy, Xexy, Yexy, XYexy,alpha)
    
    print 'done b'
    
    del A_global, d
    

    
    rhs =  np.zeros((4,n,n))
   
    b1 =  np.reshape(b[0:int_n**2], (int_n, int_n))
#    bb1 =  np.reshape(bb[0:int_n**2], (int_n, int_n))
    
#    print bb1-b1
    b2 =  np.reshape(b[int_n**2:2*int_n**2], (int_n, int_n))
    b3 =  np.reshape(b[2*int_n**2: 3*int_n**2], (int_n, int_n))
    b4 =  np.reshape(b[3*int_n**2:4*int_n**2], (int_n, int_n))
    
    rhs[0][1:-1,1:-1] = b1
    rhs[1][1:-1,1:-1] = b2
    rhs[2][1:-1,1:-1] = b3
    rhs[3][1:-1,1:-1] = b4
    
    
    u=np.zeros((4,n,n))
    
    u[0] = sin_soln(x1,y1)
    
    u[1] = sin_soln(x1,y1)
    
    u[2] = sin_soln(x1,y1)
    
    u[3] = sin_soln(x1,y1)
    
    
    s1 = 10
    s2 = 1
    
    rnorm=[np.linalg.norm(residue_sqrt_alpha(rhs, u, alpha, Ama)[0,1:-1,1:-1])*h] 
    enorm = [np.linalg.norm((u[0]-c(x1, y1))[1:-1,1:-1])*h]
    
    smoother_plot(c(x1, y1))
    
    for cycle in range(1, num_cyc+1):
        
#        print rhs ,'rhs'
        
#        print u, 'begin'

        u = VCycle(u,rhs, s1, s2, alpha, num_data)
#        print residue_sqrt_alpha(rhs, u, alpha, Ama)[0], rhs[0]
        rnorm.append(np.linalg.norm(residue_sqrt_alpha(rhs, u, alpha, Ama)[0,1:-1,1:-1])*h) 
        enorm.append(np.linalg.norm((u[0]-c(x1,y1))[1:-1,1:-1])*h)
        
    convergence_plot(num_cyc,rnorm)
    print u[0]
    smoother_plot(u[0])
        
    
    
    
    
MG(5,1e-8,100,2)