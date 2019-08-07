#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 10:48:16 2018

@author: shilu
"""

import numpy as np
from numpy import  pi, sin, cos, exp, inf
from scipy.sparse.linalg import spsolve
from scipy import eye, zeros, linalg
from numpy import linalg as LA
from copy import copy

from MGsolverUp import MGsolverUP

from grid.Grid import Grid

from BuildSquare import build_square_grid
from BuildEquation import build_equation_linear_2D, set_polynomial_linear_2D,\
Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY, NIn_triangle, build_matrix_fem_2D


from grid.function.FunctionStore import zero, exp_soln, exp_rhs, sin_soln, sin_rhs, linear, plain, l2, l3, sin4


def JacobiDown(u, rhs, totalrhs, levelnumber, dataX, dataY, data, alpha):
    
    """ Block Jacobi Method. On each level of grid (same size as initial grid), invoke corresponding matrices
    
    A, L, G1, G2, d and boundaries h1, h2, h3, h4 
    
    """
    
    print u[0], 'jacobi'
    from TPSFEM import Amatrix, Lmatrix, G1, G2, dvector
    
    # Set entire FEM grid on the level 
    grid = Grid()
    
    #print levelnumber
    # The size of FEM grid that should be in alignment with the size grid on the each level
    n =2 ** (levelnumber+1)+1

    h=1/float(n-1)
    
    # Set the mesh grid, that is interior
    #x1, y1 = np.meshgrid(np.arange(h, 1, h), np.arange(h, 1, h))
    x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    
    
    
    # Calculating the spacing
    #h=1/float(n-1)
    
    # 
    #x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    
    true_soln = zero


    
    # Build the square
    build_square_grid(n, grid, true_soln)
    
    # Store the matrices on grid
    build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  dataX, dataY)
    

    # Import matrices on each level of grid
    Amatrix = Amatrix(grid)

    
    Lmatrix =Lmatrix(grid).todense()
    #print Lmatrix
    

    
    
    dvector = dvector(grid, data)

    
    G1 = G1(grid)
    
    G2 = G2(grid)
    

    
    # Get the current size of RHS function 
    [xdim,ydim] = rhs[0][1:-1, 1:-1].shape
  

    c=np.reshape(u[0][1:-1,1:-1], (xdim**2,1))

    g1 = np.reshape(u[1][1:-1,1:-1], (xdim**2,1))

    
    g2 = np.reshape(u[2][1:-1,1:-1], (xdim**2,1))
    #print g2, 'g2'
    
    w = np.reshape(u[3][1:-1,1:-1], (xdim**2,1))
    #print c, 'CCCCCCC'
 
##############################################################################
# Sparse solver
    
     #Solving LC = f1+ G1g1+ G2g2 to get new c-grid 
    crhs = zeros([xdim+2, ydim+2])
    
    #print rhs[0].shape[0], G1.shape[0], g1.shape[0]
    
    crhs[1:-1, 1:-1] = (rhs[0][1:-1,1:-1]+ np.reshape(G1 * g1+ G2 * g2, (xdim, ydim)))

    totalrhs[0]=copy(crhs) 


    u[0][1:-1, 1:-1] = np.reshape(spsolve(Lmatrix, np.reshape(crhs[1:-1,1:-1],(xdim**2, 1))), (xdim, ydim))
    
    
    
    # Get new g1-grid

    g1rhs = zeros([xdim+2, ydim+2])

    g1rhs[1:-1, 1:-1] = (rhs[1][1:-1,1:-1]+ np.reshape( G1.T *  w ,(xdim, ydim)))

    totalrhs[1]=copy(g1rhs) / float(alpha)

    u[1][1:-1, 1:-1] = np.reshape(spsolve(alpha*Lmatrix, np.reshape(g1rhs[1:-1,1:-1],(xdim**2, 1))), (xdim, ydim))



    
    # Get new g2-grid

    g2rhs = zeros([xdim+2, ydim+2])

    g2rhs[1:-1,1:-1] = (rhs[2][1:-1,1:-1]+ np.reshape( G2.T *  w ,(xdim, ydim)))

    totalrhs[2]=copy(g2rhs) / float(alpha)
    
    u[2][1:-1, 1:-1] = np.reshape(spsolve(alpha*Lmatrix, np.reshape(g2rhs[1:-1,1:-1],(xdim**2, 1))), (xdim, ydim))
    

    
    # Get new w-grid

    wrhs = zeros([xdim+2, ydim+2])

    wrhs[1:-1,1:-1] = (rhs[3][1:-1,1:-1]- np.reshape(Amatrix * c, (xdim, ydim)))

    totalrhs[3]=copy(wrhs)
    u[3][1:-1, 1:-1] = np.reshape(spsolve(Lmatrix, np.reshape(wrhs[1:-1,1:-1],(xdim**2, 1))), (xdim, ydim))

    
    
    
    
###############################################################################
# Multigrid solver
    
    
#   
#    # Initialise rhs for c
#    crhs = zeros([xdim+2, ydim+2])
#    
#    # Compute  interior crhs (f1 + G1*g1 + G2*g2) * h**2
#    crhs[1:-1, 1:-1] = (rhs[0][1:-1,1:-1]+ np.reshape(G1 * g1+ G2 * g2, (xdim, ydim)))*float(h**2)
#    
#    # Store new crhs
#    totalrhs[0]=copy(crhs)/float(h**2)
#    
#    # Solve for c
#    u[0] = MGsolverUP(5, crhs, u[0], levelnumber)
#    
#    # Intialise rhs for new g1
#    g1rhs = zeros([xdim+2, ydim+2])
#    
#    # Compute  interior g1rhs (f2 + G1.T * w)  * h**2 / alpha
#    g1rhs[1:-1, 1:-1] = (rhs[1][1:-1,1:-1]+ np.reshape( G1.T *  w ,(xdim, ydim)))*float(h**2)/float(alpha)
#    
#    # Store g1rhs 
#    totalrhs[1]=copy(g1rhs)/float(h**2)
#    
#    # Solve for new g1
#    u[1] = MGsolverUP(5, g1rhs, u[1], levelnumber)
#    
#    # Initialise rhs for g2
#    g2rhs = zeros([xdim+2, ydim+2])
#    
#    # Compute interior g2rhs (f3 + G2.T * w)** h**2 / alpha
#    g2rhs[1:-1,1:-1] = ((rhs[2][1:-1,1:-1]+ np.reshape( G2.T *  w ,(xdim, ydim))))*float(h**2)/float(alpha)
#    
#    # Store g2rhs
#    totalrhs[2]=copy(g2rhs)/float(h**2)
#    
#    # Solve for new g2
#    u[2] = MGsolverUP(5, g2rhs, u[2], levelnumber)
#    
#    # Initialise rhs for w
#    wrhs = zeros([xdim+2, ydim+2])
#    
#    # Compute interior wrhs (f4 - Ac) * h**2
#    wrhs[1:-1,1:-1] = (rhs[3][1:-1,1:-1]- np.reshape(Amatrix * c, (xdim, ydim)))*float(h**2)
#    
#    # Store wrhs
#    totalrhs[3]=copy(wrhs)/float(h**2)
#    
#    # Solve for new w
#    u[3] = MGsolverUP(5, wrhs, u[3], levelnumber)
#    
#    
#    print totalrhs[0], 'TOTALRIGHT'
#    print u[0], 'g1after'
#    
    return u
    