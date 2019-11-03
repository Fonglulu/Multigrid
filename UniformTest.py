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
Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY, In_triangle, build_matrix_fem_2D

from Triangle import triangle_iterator, Interior_triangle_iterator
from grid.NodeTable import node_iterator, not_slave_node
from grid.ConnectTable import connect_iterator 

from BuildSquare import build_square_grid_matrix
from grid.function.FunctionStore import zero, exp_soln, exp_rhs, sin_soln, sin_rhs, linear, plain, l2, l3, sin4, exy, cos2


from mpl_toolkits.mplot3d import Axes3D

from PlotPackman import plot_fem_grid
from operator import itemgetter
from copy import copy
from functions import Linear, Xlinear, Ylinear, Linear2x, Xlinear2x, Ylinear2x, Zero, Exy, Xexy, Yexy, XYexy, Cube, Xcube, Ycube, XYcube
from Test_Matrices  import Amatrix,Stencil_A, Lmatrix, Stencil_L,  G1, G2, FD_G1, FD_G2, h1_bd, h2_bd, h3_bd, h4_bd, dv



def Setup_Grid(i):
    #from quick_d import dvector
    
    global grid, Coord, Nl,  h, Crhs, g1rhs, g2rhs, wrhs, data, intx, inty,nodes,quick_d


    x = np.linspace(0, 1.0,10)
    y = np.linspace(0, 1.0,10)
    X, Y = np.meshgrid(x,y)
    
    Z2 =Exy(X,Y)
    
    
    data = Z2.flatten()
    coordx = X.flatten()
    coordy = Y.flatten()
    Coord = zip(coordx, coordy)
    
    
    
    grid = Grid()
    
    n = 2**i+1
    
    h =1/float(n-1)
    
    intx, inty = np.meshgrid(np.arange(h, 1, h), np.arange(h, 1, h))

    x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    nodes = np.vstack([x1.ravel(), y1.ravel()]).T
    
#    quick_d = dvector(Coord, data, nodes, n)/float(len(Coord))
#    
#    quick_d = np.reshape(quick_d, (n,n))[1:-1,1:-1]
#    
#    quick_d = np.reshape(quick_d, ((n-2)**2, 1))
    
    true_soln = zero
    
    # Boundares
    Crhs = Exy

    g1rhs = Xexy

    g2rhs = Yexy

    wrhs = XYexy

    # Build node-eged structure that allows the connection 
    build_square_grid(n, grid, true_soln) 
    
    
    # Store the matrices values on grid nodes
    build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)
    
    # Set value for only interior node.
    
    Nl=[]
    
    for node in (not_slave_node(grid)):
        
        
    
        Nl.append([node,node.get_node_id()])
    
    #print 'Start sort'
    #Sort by node id.
    Nl=sorted(Nl, key = itemgetter(1))
    #print 'Finish sort'
    
    #print 'START MAKE_LIST'

    Nl=[node[0] for node in Nl]
    
    #print 'END MAKE_LIST'
    
    #print 'START SETVALUE'
        
    for node in Nl:
        node.set_value(Nl.index(node))
        #print node._coord, node.get_value()
        
    #print 'END SETVALUE'
    
    
#    Nl = 0
#    ms of the time and effort Sam and me putting into making workshop materials. We decide to talk to admin to have additional wages. Can you support us for the case
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
    


    
    Ama = Amatrix(grid, Coord, Nl, h).todense()
#   
#    
#    Lma = Lmatrix(grid, Nl).todense()
    
    
    #Ama = Stencil_A(grid, Nl,h).todense()

    

    Lma = Stencil_L(grid, Nl).todense()

    
    if discretisation == 'FEM':
        

    
        G1ma = G1(grid, Nl)#.todense()

        G2ma = G2(grid, Nl)#.todense()

        
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
    
    S_constraint = bmat([[np.diag(np.diag(Ama)), None, None,  math.sqrt(alpha)*Lma],\
                        [None, np.diag(np.diag(Lma)), None, -G1ma.T],\
                        [None, None, np.diag(np.diag(Lma)),  -G2ma.T],\
                        [ math.sqrt(alpha)*Lma, -G1ma, -G2ma, None]])
    
    
    
    S_diag = bmat([[math.sqrt(alpha)*Lma, -G1ma, -G2ma, None],\
                       [None, Lma, None, -G1ma.T],\
                       [None, None, Lma,  -G2ma.T],\
                       [Ama, None, None, math.sqrt(alpha)*Lma]])
    

    
    S = bmat([[Ama, None, None,  None],\
                        [None, Lma, None, -G1ma.T],\
                        [None, None, Lma,  -G2ma.T],\
                        [ None, -G1ma, -G2ma, None]])
    
    Stoke = bmat([[Lma, None, -G1ma.T],\
                        [ None, Lma,  -G2ma.T],\
                        [  -G1ma, -G2ma, None]])
    
    Test = bmat([[Ama, None, None,  np.identity(Ama.shape[0])],\
                        [None, Lma, None, -G1ma.T],\
                        [None, None, Lma,  -G2ma.T],\
                        [  np.identity(Ama.shape[0]), -G1ma, -G2ma, None]])
    
    B = bmat([[math.sqrt(alpha)*Lma, -G1ma, -G2ma, None]])
    

    
    



    
    return   S_sym, S.todense(), Stoke.todense(), Test.todense(), Ama


def Setup_Rhs(alpha,i):
    
    global grid, Coord, Nl, Chrs, g1rhs, g2rhs, wrhs,quick_d
    
    Setup_Grid(i)

   
    h1 = h1_bd(grid, Crhs, wrhs, alpha, Nl, Coord)
    h2 = h2_bd(grid, g1rhs, wrhs, alpha, Nl)
    h3 = h3_bd(grid, g2rhs, wrhs, alpha, Nl)
    h4 = h4_bd(grid, Crhs, g1rhs, g2rhs, Nl)
    d  = dv(grid, data, Nl, Coord)
    
    
    rhs = zeros([4*len(Nl),1])
    rhs[0:len(Nl)] = d - h1
    rhs[len(Nl):2*len(Nl)] = -h2/float(math.sqrt(alpha))
    rhs[2*len(Nl):3*len(Nl)] = -h3/float(math.sqrt(alpha))
    rhs[3*len(Nl):] = -h4*float(math.sqrt(alpha))
    
    
    rhs_test = zeros([4*len(Nl),1])
    rhs_test[0:len(Nl)] =  d
    rhs_test[len(Nl):2*len(Nl)] = -h2
    rhs_test[2*len(Nl):3*len(Nl)] = -h3
    rhs_test[3*len(Nl):] = -h4
    rhs_test = np.reshape(rhs_test,(4,int(math.sqrt(len(Nl))), int(math.sqrt(len(Nl)))))
    
    
    
    return rhs,rhs_test,h3
    
    


def Ploteig(precond_matrix, matrix, section):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import subplots
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from pylab import title
    from numpy.linalg import inv
    from copy import copy
    
    global Nl
    
    #Gram = np.matrix(matrix.T)*np.matrix(matrix)
##    Gram2 =  np.matmul(matrix.T, matrix)
##    Gram3 = Gram-Gram2
#                        
#    Diag =np.diag(Gram)  # Take the diagonal 
#    Diag =copy(Diag)
#    
##    d1 = np.array([1 for i in Diag[0:1*len(Nl)]])
###    
###    Diag[0*len(Nl):1*len(Nl)]=d1
##    diag = np.asarray([k for k in Diag])
##    
##    #
#    Diagonal = inv(np.matrix(np.diag(Diag)))
##    #
#    Diagonal = np.matrix(np.sqrt(Diagonal))
##    print LA.cond(Diagonal)
##    
##    #
#    precond_Gram = Diagonal* Gram *Diagonal
    
    
    matrix = np.matmul(np.matrix(inv(precond_matrix)), np.matrix(matrix))
    
    #w eigenvalues, v eigenvectors
    Eigenvalue, Eigenvector=LA.eig(matrix)
    #eigvector, eigvalue, vh=np.linalg.svd(Amatrix)
    print 'max', max(abs(Eigenvalue)), 'min', min(abs(Eigenvalue)), 'cond', LA.cond(matrix)
    
    abs_eigvalue=abs(Eigenvalue)
    abs_eigvalue=abs_eigvalue.tolist()
    
    number = number_of_one(abs_eigvalue)
    
    print number

    fig1 = plt.figure(figsize=(10,10))
    ax1 = fig1.add_subplot(2,1,1)
    ax1.scatter(range(len(abs_eigvalue)),sorted(abs_eigvalue), s=1)
    
    
    ax2 = fig1.add_subplot(2,1,2)
    ax2.scatter(Eigenvalue.real, Eigenvalue.imag,s=1)
    #title('eigenvalue distribution')
    plt.show()
    
    
 

    
    # Return the position index. reverse gives order from big to small
    Maxeigvalue = abs_eigvalue.index(sorted(abs_eigvalue, reverse= True)[-2])
    
    # Find the eigenvector corresponding to max eigenvalue
    Maxeigvector = abs(Eigenvector[:,Maxeigvalue])

    Maxeigvector = Maxeigvector.tolist()
    
    # eigenvector plot in 2D
    plt.scatter(range(len(Maxeigvector)),Maxeigvector, s=1)
    title('eigenvector plot of assigned eigenvalue')
    plt.show()
    
    
    if section == 'c':
        
        
        s = 0
        
    elif section == 'g1':
        
        s = 1
        
    elif section == 'g2':
        
        s=2
        
    elif section == 'w':
        
        s=3
        
    
    # eigenvector plot in 2D for a section
    plt.scatter(range(len(Maxeigvector[s*len(Nl):(s+1)*len(Nl)])),Maxeigvector[s*len(Nl):(s+1)*len(Nl)], s=1)
    title('Section eigenvector plot of assigned eigenvalue')
    plt.show()
                    
    #Take a section of eigenvector for 
    Maxeigvector= Maxeigvector[s*len(Nl):(s+1)*len(Nl)]
    print Maxeigvector[0:10]
                        
                    
                        
    # eigenvector plot in 3D
    Maxeigvector = np.reshape(Maxeigvector, (int(np.sqrt(len(Nl))), int(np.sqrt(len(Nl)))))
                        
    ax = Axes3D(plt.gcf())
                            
    ax.plot_surface(intx,inty, Maxeigvector)  
    
    title('Section eigenvector plot of assigned eigenvalue')
            






def number_of_one(v):
    
    count = 0
    
    for i in v:
        
        if abs(i-1) < 1e-8:
            
            count +=1
    return count
            


