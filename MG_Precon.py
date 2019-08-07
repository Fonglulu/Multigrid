#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 11:30:45 2019

@author: shilu
"""

###############################################################################
# This routine defines multigrid preconditioning step as LinearOperator
###############################################################################

import numpy as np
from scipy.sparse.linalg import LinearOperator
#from QMR import alpha

alpha = 1
def residue_sqrt_alpha(rhs, u, alpha):
    
    """ This routine takes the approximation  of u, computes the residue by
    r=Au-f, where A is the modified version with square root of alpha.
    """
    
    import numpy as np
    from Rich_uniform import Astencil, G1stencil, G2stencil, Lstencil

    # Get the current size of RHS function
    [xdim,ydim] = rhs[0][1:-1,1:-1].shape
    
    h = 1/ float(xdim+2-1)

    
    # Initialise the residual

    r=np.zeros((4, rhs.shape[1],rhs.shape[2]))
    
    r[0] = rhs[0] - np.sqrt(alpha)*Lstencil(u[0])+ G1stencil(u[1], h) + G2stencil(u[2] ,h)
    
    r[1] = rhs[1] - Lstencil(u[1]) - G1stencil(u[3],h)
    
    r[2] = rhs[2] - Lstencil(u[2]) - G2stencil(u[3], h)
    
    r[3] = rhs[3] - Astencil(u[0], h) - np.sqrt(alpha)* Lstencil(u[3])
    

    
    
    return r

def Interpolation(uc):
    
    """ 
    Interpolate coarse grid to fine grid
    
    Input: current approximation on coarse grid uc

    Output: Interpolated approximation on find grid    
    """
    [depth, xdim, ydim] = uc.shape
    
    # Initialise a next fine grid
    xnodes = 2*xdim-1
    ynodes = 2*ydim-1
    uf = np.zeros((depth, xnodes,ynodes))
    
    
    # For even ordered i and j on fine grid
    for k in range(depth):
        for i in range(xdim):
            for j in range (ydim):
                uf[k, 2*i, 2*j]=uc[k, i,j]
    

    # For even ordered j on fine grid on fine grid
    for k in range(depth):
        for i in range(0, ynodes, 2):
            for j in range(1, xnodes-1, 2):
                uf[k,i,j]=0.5*(uf[k,i,j-1]+uf[k,i,j+1])

        
    # For even ordered i on fine grid on fine grid
    for k in range(depth):
        for i in range(1, xnodes-1, 2):
            for j in range (0, ynodes, 2):
                uf[k,i,j]=0.5*(uf[k,i-1,j]+uf[k,i+1,j])
    
    # For odd ordered i and j on fine grid on fine grid
    for k in range(depth):
        for i in range (1, xnodes-1, 2):
            for j in range (1, ynodes-1, 2):
                uf[k,i,j]=0.25*(uf[k,i-1,j]+uf[k,i+1,j]+uf[k,i,j-1]+uf[k,i,j+1])#    

            
            
    return uf




def FW_Restriction(uf):
    
    """ 
    Restrict fine grid to coarse grid by full weighting operator
    
    Input: current approximation on fine grid uf
    
    Output: Restricted approximation on coarse grid 
    
    """

    # Get the current fine grid
    [depth, xdim, ydim] = uf.shape

    
    # Coarse grid size
    xnodes = int((xdim+1)/2)
    ynodes = int((ydim+1)/2)
    
    # Set coarse grid
    uc = np.zeros((depth, xnodes,ynodes))
    
    
    # Find the values from the original positions
    for k in range(depth):
        for i in range(1, xnodes-1):
            for j in range(1, ynodes-1):
                
                uc[k,i,j] = uf[k,2*i,2*j] + 0.5*(uf[k, 2*i-1, 2*j]+uf[k, 2*i+1, 2*j] + uf[k,2*i, 2*j-1]+\
                    uf[k,2*i,2*j+1]) + 0.25* (uf[k,2*i-1,2*j-1] + uf[k,2*i-1, 2*j+1] + uf[k,2*i+1, 2*j-1] + uf[k,2*i+1, 2*j+1])
                
    return uc




def VCycle(u, v,  s1, s2, alpha):
    """ This routine implements the recurssion version of V-cycle
        input: current approximation of u
               rhs 
               s1: number of iterations for relaxtions applied downwards
               s2: number of iterations for relaxtions applied upwards
               alpha
        outpu: new approximation of u
    """
    
    from Rich_uniform import Rich
    
    
    if u[0].shape[0] != 3:
        
        for sweeps in range(s1):
            #print 'downwards'
            u = Rich(u, v, alpha)
        
        new_rhs = residue_sqrt_alpha(v, u, alpha)

        
        new_rhs = FW_Restriction(new_rhs)
        
        uc = np.zeros((4, new_rhs[0].shape[0], new_rhs[0].shape[0]))
        
        uc = VCycle(uc, new_rhs, s1, s2, alpha)
        
        u = u + Interpolation(uc)
        
        
    
    for sweeps in range(s2):
        
             u = Rich(u, v, alpha)
        
    return u


def multigrid_matvec(rhs):
    
    """ Compute A^-1 v by solving Ax=v with multigrid method"""
    
    rhs_list = rhs.tolist()
    
    # size of interior mesh (no boundary)
    n = np.int(np.sqrt(len(rhs_list)/4))
    
    h1 = np.reshape(rhs_list[0:n**2], (n,n))
    
    h2 = np.reshape(rhs_list[n**2:2*n**2], (n,n))
    
    h3 = np.reshape(rhs_list[2*n**2:3*n**2], (n,n))
    
    h4 = np.reshape(rhs_list[3*n**2:4*n**2], (n,n))
    
#    rhs = np.zeros((4,n,n))
#    
#    rhs[0] = h1
#    
#    rhs[1] = h2
#    
#    rhs[2] = h3
#    
#    rhs[3] = h4
    
    rhs=np.zeros((4,n+2,n+2))
    
    rhs[0, 1:-1, 1:-1] = h1
    
    rhs[1, 1:-1, 1:-1] = h2
    
    rhs[2, 1:-1, 1:-1] = h3
    
    rhs[3, 1:-1, 1:-1] = h4
    
    
    
#    layer,n, n = np.shape(rhs) 
    
    u=np.zeros((4,n+2,n+2))
    
    #u=np.zeros((4,n,n))
    
    s1 = 10
    s2 = 10
    
    for cycle in range(1,10):
        
         u = VCycle(u,rhs, s1, s2, alpha)
         
         
    u = u[:,1:-1,1:-1]
         
    u = u.reshape((4*(n)**2,))
    #print np.shape(u)
    return u



def pre_multigrid(n):
    
    M = LinearOperator((4*n**2,4*n**2), matvec= multigrid_matvec, rmatvec=multigrid_matvec)
    print 4*n**2
    return M
    