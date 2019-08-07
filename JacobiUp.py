#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 10:54:28 2018

@author: shilu
"""

import numpy as np
from numpy import  pi, sin, cos, exp, inf
from scipy.sparse.linalg import spsolve
from scipy import eye, zeros, linalg
from numpy import linalg as LA
from copy import copy


from grid.Grid import Grid
from MGsolverUp import MGsolverUP


from BuildSquare import build_square_grid
from BuildEquation import build_equation_linear_2D, set_polynomial_linear_2D,\
Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY, NIn_triangle, build_matrix_fem_2D

from grid.function.FunctionStore import zero, linear


def JacobiUp(u, rhs, totalrhs, levelnumber, dataX, dataY, data, alpha):
    
    """ Block Jacobi Method. On each level of grid (same size as initial grid), invoke corresponding matrices
    
    A, L, G1, G2, d and boundaries h1, h2, h3, h4 
    
    """
    
    print u[0], 'Upjacobi'
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
    #print Amatrix

    
    
    dvector = dvector(grid, data)
    #print dvector
    
    #Lmatrix = Lmatrix(grid)
    #print Lmatrix.todense()
    
    
    G1 = G1(grid)

    
    G2 = G2(grid)
    

    
    # Get the current size of RHS function 
    [xdim,ydim] = rhs[0][1:-1, 1:-1].shape
    c=np.reshape(u[0][1:-1,1:-1], (xdim**2,1))
    #print c, 'c'
    
    
    #print v[1]
    g1 = np.reshape(u[1][1:-1,1:-1], (xdim**2,1))
    #print g1
    #print np.size(g1, 0), 'g1'
    
    g2 = np.reshape(u[2][1:-1,1:-1], (xdim**2,1))
    #print g2, 'g2'
    
    w = np.reshape(u[3][1:-1,1:-1], (xdim**2,1))
    #print w, 'w'
    
    
    # Solving LC = f1+ G1g1+ G2g2 to get new c-grid 
    #crhs = zeros([xdim+2, ydim+2])
    crhs = 300*np.ones((xdim+2,ydim+2))
    #print crhs.shape[0]
    #rhs[0][1:-1,1:-1] = rhs[0][1:-1,1:-1]+ np.reshape(G1 * g1+ G2 * g2, (xdim, ydim))

    crhs[1:-1, 1:-1] = (rhs[0][1:-1,1:-1]+ np.reshape(G1 * g1+ G2 * g2, (xdim, ydim)))*float(h**2)
    #print totalrhs[0].shape[0], 'TOTALDIM'
    totalrhs[0]=copy(crhs) 
    #crhs = np.reshape(rhs[0][1:-1,1:-1], (xdim**2,1))+ G1 * g1+ G2 * g2
    #print crhs, 'crhs'
    
    
    # Reshape approximation of c so that corresponds to nodes
    
    
    # Applying single multigrid solver
    #print crhs, 'CRHSSSSS'
    u[0] = MGsolverUP(5, crhs, u[0], levelnumber)
    #u[0][1:-1,1:-1] = np.reshape(spsolve(Lmatrix, crhs), (xdim, ydim))
    
    
    
    # Get new g1-grid
   # g1rhs = (np.reshape(rhs[1][1:-1,1:-1], (xdim**2,1)) + (G1.T)* w)/float(alpha)
    #g1rhs = zeros([xdim+2, ydim+2])
    g1rhs = 300*np.ones((xdim+2,ydim+2))
    g1rhs[1:-1, 1:-1] = ((rhs[1][1:-1,1:-1]+ np.reshape( G1.T *  w ,(xdim, ydim)))*float(h**2))/float(alpha)
    #rhs[1][1:-1,1:-1] = (rhs[1][1:-1,1:-1]+ np.reshape( G1.T *  w ,(xdim, ydim)))/float(alpha)
    totalrhs[1]=copy(g1rhs) 
    u[1] = MGsolverUP(5, g1rhs, u[1], levelnumber)
    #print u[1], 'g1after'
    #u[1][1:-1,1:-1] = np.reshape(spsolve(Lmatrix, g1rhs), (xdim, ydim))
    #print u[1], 'u[1]'

    #print rhs[1], 'rhs[1]'
    
    # Get new g2-grid
    #g2rhs = (np.reshape(rhs[2][1:-1,1:-1], (xdim**2,1))+ (G2.T) * w)/float(alpha)
    #g2rhs = zeros([xdim+2, ydim+2])
    g2rhs = 300*np.ones((xdim+2,ydim+2))
    #g2rhs[1:-1,1:-1] = ((rhs[2][1:-1,1:-1]+ np.reshape( G2.T *  w ,(xdim, ydim)))*float(h**2))/float(alpha)
    g2rhs[1:-1,1:-1] = (rhs[2][1:-1,1:-1]+ np.reshape( G2.T *  w ,(xdim, ydim)))/float(alpha)
    #rhs[2][1:-1,1:-1] = (rhs[2][1:-1,1:-1]+ np.reshape( G1.T *  w ,(xdim, ydim)))/float(alpha)
    #u[2][1:-1,1:-1] = np.reshape(spsolve(Lmatrix, g2rhs), (xdim, ydim))
    totalrhs[2]=copy(g2rhs)
    u[2] = MGsolverUP(5, g2rhs, u[2], levelnumber)
    
    # Get new w-grid
#    wrhs = dvector - h1 - np.matmul(Amatrix, c)
#    wrhs = np.reshape(wrhs, (xdim,ydim))
    #print rhs[3][1:-1,1:-1]
    #wrhs = zeros([xdim+2, ydim+2])
    wrhs = 300*np.ones((xdim+2,ydim+2))
    wrhs[1:-1,1:-1] = (rhs[3][1:-1,1:-1]- np.reshape(Amatrix * c, (xdim, ydim)))#*float(h**2)
    #rhs[3][1:-1,1:-1] = rhs[3][1:-1,1:-1]- np.reshape(Amatrix * c, (xdim, ydim))
    #wrhs = np.reshape(rhs[3][1:-1,1:-1], (xdim**2,1))- Amatrix * c
    totalrhs[3]=copy(wrhs) 
    
    u[3] = MGsolverUP(5, wrhs, u[3], levelnumber)
    #u[3][1:-1,1:-1] = np.reshape(spsolve(Lmatrix, np.reshape(wrhs, (xdim**2,1)), (xdim, ydim)))
    
    #print u[0], 'g1after'
    
    return u
    