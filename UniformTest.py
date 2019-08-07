#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:11:17 2019

@author: shilu
"""

import timeit
import numpy as np
from numpy import  pi, sin, cos, exp, inf
from scipy import zeros, linalg
from scipy.sparse import csc_matrix, lil_matrix, bmat, coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from numpy import linalg as LA
import math



from grid.Grid import Grid
from BuildSquare import build_square_grid
from BuildEquation import build_equation_linear_2D, set_polynomial_linear_2D,\
Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY, NIn_triangle, build_matrix_fem_2D

from Triangle import triangle_iterator, Interior_triangle_iterator
from grid.NodeTable import node_iterator, not_slave_node
from grid.ConnectTable import connect_iterator 

from BuildSquare import build_square_grid_matrix
from grid.function.FunctionStore import zero, exp_soln, exp_rhs, sin_soln, sin_rhs, linear, plain, l2, l3, sin4, exy, cos2
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

from PlotPackman import plot_fem_grid
from operator import itemgetter
from copy import copy
from functions import Linear, Xlinear, Ylinear, Zero, Exy, Xexy, Yexy, XYexy, Cube, Xcube, Ycube, XYcube
from Test_Matrices  import Amatrix,Stencil_A, Lmatrix, Stencil_L,  G1, G2, FD_G1, FD_G2, h1_bd, h2_bd, h3_bd, h4_bd, dvector



def Setup_Grid(i):
    
    global grid, Coord, Nl,  h, Crhs, g1rhs, g2rhs, wrhs, data, x1, y1


    x = np.linspace(0, 1.0,100)
    y = np.linspace(0, 1.0,100)
    X, Y = np.meshgrid(x,y)
    
    Z2 = Linear(X,Y)
    
    
    data = Z2.flatten()
    coordx = X.flatten()
    coordy = Y.flatten()
    Coord = zip(coordx, coordy)
    
    
    
    grid = Grid()
    
    n = 2**i+1
    
    h =1/float(n-1)
    
    x1, y1 = np.meshgrid(np.arange(0, 1-h, h), np.arange(0, 1-h, h))
    
    
    true_soln = zero
    
    # Boundares
    Crhs = Linear

    g1rhs = Xlinear

    g2rhs = Ylinear

    wrhs = Zero


    build_square_grid(n, grid, true_soln)
    
    build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)
    
    # Set value for only interior node.
    
    Nl=[]
    
    for node in (not_slave_node(grid)):
    
        Nl.append([node,node.get_node_id()])
    
    print 'Start sort'
    Nl=sorted(Nl, key = itemgetter(1))
    print 'Finish sort'
    
    print 'START MAKE_LIST'

    Nl=[node[0] for node in Nl]
    
    print 'END MAKE_LIST'
    
    print 'START SETVALUE'
        
    for node in Nl:
        node.set_value(Nl.index(node))
        
    print 'END SETVALUE'
    
    
#    Nl = 0
#    
#    for node in node_iterator(grid):
#        
#        if not node.get_slave():
#            
#            node.set_value(Nl)
#            
#            Nl = Nl+1
    
    

def Setup_Matrices(alpha, i, discretisation):
    
    global grid, Coord, Nl,  h
    
    Setup_Grid(i)
    


    
#    Ama = Amatrix(grid, Coord, Nl, h).todense()
#   
#    
#    Lma = Lmatrix(grid, Nl).todense()
    
    
    Ama = Stencil_A(grid, Nl,h).todense()

    

    Lma = Stencil_L(grid, Nl).todense()

    
    if discretisation == 'FEM':
        

    
        G1ma = G1(grid, Nl)

        G2ma = G2(grid, Nl)

        
    else:
        
        print 'G1-START'
        G1ma = FD_G1(grid, Nl, h).todense()
        print 'G1-FINISH'
        print 'G2-START'
        G2ma = FD_G2(grid, Nl, h).todense()
        print 'G2-FINISH'
        
    #ZeroMatrix = csr_matrix((len(Nl),len(Nl)))
    
    #G1ma = ZeroMatrix
    #G2ma = ZeroMatrix
    
    
    S_sym = bmat([[Ama, None, None,  math.sqrt(alpha)*Lma],\
                        [None, Lma, None, -G1ma.T],\
                        [None, None, Lma,  -G2ma.T],\
                        [ math.sqrt(alpha)*Lma, -G1ma, -G2ma, None]])
    
    
    
    S_diag = bmat([[math.sqrt(alpha)*Lma, -G1ma, -G2ma, None],\
                       [None, Lma, None, -G1ma.T],\
                       [None, None, Lma,  -G2ma.T],\
                       [Ama, None, None, math.sqrt(alpha)*Lma]])
    
    



    
    return   S_sym.todense(), S_diag.todense() # SymBigMat.todense()  #SymBigMat SymBigMat.todense()


def Setup_Rhs(alpha,i):
    
    global grid, Coord, Nl, Chrs, g1rhs, g2rhs, wrhs
    
    Setup_Grid(i)

    
    h1 = h1_bd(grid, Crhs, wrhs, Nl, Coord)
    h2 = h2_bd(grid, g1rhs, wrhs, alpha, Nl)
    h3 = h3_bd(grid, g1rhs, wrhs, alpha, Nl)
    h4 = h4_bd(grid, Crhs, g1rhs, g2rhs, Nl)
    d  = dvector(grid, data, Nl, Coord)
    
    
    rhs = zeros([4*len(Nl),1])
    rhs[0:len(Nl)] = d - h1
    rhs[len(Nl):2*len(Nl)] = -h2/float(math.sqrt(alpha))
    rhs[2*len(Nl):3*len(Nl)] = -h3/float(math.sqrt(alpha))
    rhs[3*len(Nl):] = -h4*float(math.sqrt(alpha))
    
    
    
    
    return rhs
    
    


def Ploteig(matrix):
    import matplotlib.pyplot as pp
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from numpy.linalg import inv
    from copy import copy
    
    global Nl
    
    Gram = np.matrix(matrix.T)*np.matrix(matrix)
    #Gram2 =  np.matmul(matrix.T, matrix)
    #Gram3 = Gram-Gram2

#    Diag =np.diag(Gram)
#    Diag =copy(Diag)
#    
##    d1 = np.array([1 for i in Diag[0:1*len(Nl)]])
##    
##    Diag[0*len(Nl):1*len(Nl)]=d1
#    #diag = np.asarray([k for k in Diag])
#    
#    #
#    Diagonal = inv(np.matrix(np.diag(Diag)))
#    #
#    Diagonal = np.matrix(np.sqrt(Diagonal))
#    print LA.cond(Diagonal)
#    
#    #
#    Gram1 = Diagonal* Gram *Diagonal
    
    
    
    
    #w eigenvalues, v eigenvectors
    Eigenvalue, Eigenvector=LA.eigh(Gram)
    #eigvector, eigvalue, vh=np.linalg.svd(Amatrix)
    print 'max', max(Eigenvalue), 'min', min(abs(Eigenvalue)), 'cond', LA.cond(Gram)
    
    
    Leigvalue=abs(Eigenvalue)
    Leigvalue=Leigvalue.tolist()
    
    pp.scatter(range(len(Leigvalue)),sorted(Leigvalue))
    pp.show()
    
    
 

    
    # Return the position index. reverse gives order from big to small
    Lmaxeigvalue = Leigvalue.index(sorted(Leigvalue, reverse= True)[3])
    print sorted(Leigvalue, reverse= True)[3]
    
    # Find the eigenvector corresponding to max eigenvalue
    Lmaxeigvector = Eigenvector[:,Lmaxeigvalue]
    

    #Take the section of eigenvector for c
    Lmaxeigvector=Lmaxeigvector[0:len(Nl)]
    
    
    Lmaxeigvector = np.reshape(Lmaxeigvector, (int(np.sqrt(len(Nl))), int(np.sqrt(len(Nl)))))


    ax = Axes3D(plt.gcf())
    
    ax.plot_surface(x1,y1, Lmaxeigvector)   








