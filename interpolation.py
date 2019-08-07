#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 07:25:18 2018

@author: shilu
"""

import numpy as np

def Interpolation(uc):
    
    """ 
    Interpolate coarse grid to fine grid
    
    Input: current approximation on coarse grid uc

    Output: Interpolated approximation on find grid    
    """
    
    # Get the current size of approximation on coarse grid
    [xdim, ydim] = uc.shape
    
    # Initialise a next fine grid
    xnodes = 2*xdim+1
    ynodes = 2*ydim+1
    grid = np.zeros((xnodes,ynodes))
    
    
# For even ordered i and j

    for i in range(xdim):
        for j in range (ydim):
            grid[2*i+1,2*j+1]=uc[i,j]


# For even ordered j 
  
    for i in range(1, ynodes-1, 2):
        for j in range(2, xnodes-2, 2):
            grid[i,j]=0.5*(grid[i,j-1]+grid[i,j+1])

    
# For even ordered i   

    for i in range(2, xnodes-2, 2):
        for j in range (1, ynodes-1, 2):
            grid[i,j]=0.5*(grid[i-1,j]+grid[i+1,j])

# For odd ordered i and j


    for i in range(2, xnodes-2, 2):
        for j in range (2, ynodes-2, 2):
            grid[i,j]=0.25*(grid[i-1,j]+grid[i+1,j]+grid[i,j-1]+grid[i,j+1])
            
# For boundary


    for i in range(1, xnodes-1):
        grid[i, 0] = 2*grid[i,1]- grid[i,2]
        grid[i,-1] = 2*grid[i,-2]-grid[i, -3]


    for j in range(1, ynodes-1):
        grid[0, j] = 2* grid[1, j]- grid[2, j]
        grid[-1,j] = 2* grid[-2, j] - grid[-3, j]
        
# Get the current size of approximation on coarse grid
#    [depth, xdim, ydim] = uc.shape
#    
#    # Initialise a next fine grid
#    xnodes = 2*xdim+1
#    ynodes = 2*ydim+1
#    grid = np.zeros((depth, xnodes,ynodes))
#    
#    
#    # For even ordered i and j
#    for k in range(depth):
#        for i in range(xdim):
#            for j in range (ydim):
#                grid[k, 2*i+1,2*j+1]=uc[k, i,j]
#    
#
#    # For even ordered j 
#    for k in range(depth):    
#        for i in range(1, ynodes-1, 2):
#            for j in range(2, xnodes-2, 2):
#                grid[k, i,j]=0.5*(grid[k, i,j-1]+grid[k, i,j+1])
#
#        
#    # For even ordered i   
#    for k in range(depth):
#        for i in range(2, xnodes-2, 2):
#            for j in range (1, ynodes-1, 2):
#                grid[k,i,j]=0.5*(grid[k, i-1,j]+grid[k,i+1,j])
#    
#    # For odd ordered i and j
#    
#    for k in range(depth):
#        for i in range(2, xnodes-2, 2):
#            for j in range (2, ynodes-2, 2):
#                grid[k, i,j]=0.25*(grid[k, i-1,j]+grid[k, i+1,j]+grid[k, i,j-1]+grid[k, i,j+1])
#                
    # For boundary
    
#    for k in range(depth):
#        for i in range(1, xnodes-1):
#            grid[k, i, 0] = 2*grid[k, i,1]- grid[k,i,2]
#            grid[k, i,-1] = 2*grid[k,i,-2]-grid[k, i, -3]
#    
#    for k in range(depth):
#        for j in range(1, ynodes-2):
#            grid[k, 0, j] = 2* grid[k, 1, j]- grid[k,2, j]
#            grid[k, -1,j] = 2* grid[k, -2, j] - grid[k, -3, j]\
    
    
    
#    crhs =np.zeros((xdim, ydim))
#    
#    crhs = (rhs[0][1:-1,1:-1]+ np.reshape(G1 *g1 +G2* g2, (xdim, ydim)))
#    
#    u[0][1:-1,1:-1] = np.reshape(spsolve(Lmatrix, crhs), (xdim, ydim))
#    
#    
#    g1rhs = np.ones((xdim, ydim))
#    
#    g1rhs = (rhs[1][1:-1,1:-1]+ np.reshape(G1.T * w, (xdim, ydim)))
#    
#    u[1] = np.reshape(spsolve(alpha* Lmatrix, g1rhs), (xdim, ydim))
#    
#    
#    g2rhs = np.ones((xdim, xdim))
#    
#    g2rhs = (rhs[2][1:-1,1:-1]+ np.reshape(G2.T * w, (xdim, ydim)))
#    
#    u[2]  = np.reshape(spsolve(alpha* Lmatrix, g2rhs), (xdim, ydim))
#    
#    
#    w = np.ones((xdim, xdim))
#    
#    wrhs = (rhs[3][1:-1,1:-1]- np.reshape(Amatrix * c, (xdim, ydim)))
#    
#    u[3]  = np.reshape(spsolve(Lmatrix, g2rhs), (xdim, ydim))
            
            
    return grid