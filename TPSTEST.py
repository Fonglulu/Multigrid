#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 19:35:25 2018

@author: shilu
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 19:44:03 2018

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




alpha = 0.000000001


x = np.linspace(0, 1.0,50)
y = np.linspace(0, 1.0,50)
X, Y = np.meshgrid(x,y)

Z2 = Cube(X,Y)


data = Z2.flatten()
coordx = X.flatten()
coordy = Y.flatten()
Coord = zip(coordx, coordy)



grid = Grid()

i =3
n = 2**i+1

h =1/float(n-1)

x1, y1 = np.meshgrid(np.arange(0, 1-h, h), np.arange(0, 1-h, h))


true_soln = zero
#
# Boundares
Crhs = Cube
#
g1rhs = Xcube
#
g2rhs = Ycube
#
wrhs = XYcube
#
#
#
build_square_grid(n, grid, true_soln)
#
build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)







Nl=[]
for node in (not_slave_node(grid)):
    
    Nl.append([node,node.get_node_id()])
    
Nl=sorted(Nl, key = itemgetter(1))

Nl=[node[0] for node in Nl]
        
for node in Nl:
    node.set_value(Nl.index(node))


# Generate Amatrix
    
def Amatrix(grid):
    

#        node.set_value(Nl.index(node))
    
#        
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
            
    for node in Nl:
        node.set_value(Nl.index(node))
            
    Amatrix = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
        
        if not node.get_slave():
            
            idi = int(node.get_value())
            
            #print idi, node.get_coord()

            
            #print node.get_value(), node.get_node_id()._id_no
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
                idj = int(grid.get_value(endpt))
                
                #print idj, grid.get_coord(endpt)
                
                #print endpt._id_no, grid.get_value(endpt)
                
                aentry = grid.get_matrix_value(node.get_node_id(), endpt)[1]

                #print aentry
                
                if not grid.get_slave(endpt):
                    
                    Amatrix[idi, idj] =aentry


    return Amatrix/float(len(Coord))

#
Amatrix = Amatrix(grid).todense()
       
     
            
# Generate h1 boundary condition Ac + Lw = d - h1

def h1(grid, Crhs, wrhs): 
    

    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
            
    for node in Nl:
        node.set_value(Nl.index(node))
        
    
    h1 = zeros((len(Nl), 1))
        
    h = 1/float(math.sqrt(len(Nl))+1)
    
    for node in not_slave_node(grid):
        
        # Ignore slave (or boundary) nodes
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            #print node.get_value(), node.get_node_id()._id_no
            
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
                #j = int(grid.get_value(endpt))
                
                #print j
                
            
                if grid.get_slave(endpt):
                    
                    
                    coord = grid.get_coord(endpt)
            
                    c = Crhs(coord[0], coord[1])
                    
                    w = wrhs(coord[0], coord[1])

                                        
                    #print aentry
                    
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
            
                    aentry = grid.get_matrix_value(node.get_node_id(), endpt)[1]/float(len(Coord))
            
                    h1[i] += c* aentry + w * lentry
    return h1

                
h1 =h1(grid, Crhs, wrhs)


# Generate h2 boundary condition alphaLg1 -G1^Tw = -h2
    
def h2(grid, g1rhs, wrhs, alpha):
    
#    g1rhs = zero
#    
#    wrhs =zero
    
#    build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)
    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
            
    for node in Nl:
        node.set_value(Nl.index(node))
    
    
    h2= zeros([len(Nl), 1])
    
    for node in not_slave_node(grid):
        
        # Ignore slave (or boundary) nodes
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
            
                if grid.get_slave(endpt):
                    
                    coord = grid.get_coord(endpt)
                    
                    #print coord
            
                    g1 = g1rhs(coord[0], coord[1])
                    #print g1
                    
                    w = wrhs(coord[0], coord[1])

                    
                    #print endpt._id_no
            
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
                    
                    g1entry = grid.get_matrix_value(node.get_node_id(), endpt)[2]
                    
                    #print aentry
            
                    #h2[i] += alpha * g1 * lentry - w * g1entry #G1, G2
                    h2[i] += alpha * g1 * lentry +w * g1entry # -G1, -G2
                    
    return h2


h2 = h2(grid, g1rhs, wrhs, alpha)


# Generate h3 boundary condition alphaLg2 -G2^Tw = -h3

def h3(grid, g2rhs, wrhs, alpha):
    
#    g2rhs = zero
#    
#    wrhs = zero
#    
#    build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)
#    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
 
    for node in Nl:
        node.set_value(Nl.index(node))
    
    
    h3= zeros([len(Nl), 1])
    
    for node in not_slave_node(grid):
        
        # Ignore slave (or boundary) nodes
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
            
                if grid.get_slave(endpt):
                    
                    coord = grid.get_coord(endpt)
            
                    g2 = g2rhs(coord[0], coord[1])
                    
                    w = wrhs(coord[0], coord[1])
                    
                    #print c, w
                    
                    #print endpt._id_no
            
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
                    
                    # G2 entries on boundary
                    g2entry = grid.get_matrix_value(node.get_node_id(), endpt)[3]
                    
                    #print aentry
                    
                    #h3[i] += alpha * g2* lentry - w* g2entry # G1, G2
                    h3[i] += alpha * g2* lentry + w* g2entry  # -G1, -G2
                    
    return h3

h3 = h3(grid, g2rhs, wrhs, alpha)

# Generate h4 boundary condition for Lc -G1g1 -G2g2 = -h4 on boundary
    
def h4(grid, Crhs, g1rhs, g2rhs):
    
#    Crhs = plain
#    
#    g1rhs = zero
#    
#    g2rhs = zero
#
#    build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)
    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
    
    for node in Nl:
        node.set_value(Nl.index(node))
    
    
    h4 = zeros([len(Nl),1])
    
    for node in not_slave_node(grid):
    
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
            
                if grid.get_slave(endpt):
                    
                    coord = grid.get_coord(endpt)
                    
            
                    c = Crhs(coord[0], coord[1])
                    
                    g1 = g1rhs(coord[0], coord[1])
                    
                    g2 = g2rhs(coord[0], coord[1])
                    
                    #print c
                    
                    #print endpt._id_no
            
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
                    
                    g1entry = grid.get_matrix_value(node.get_node_id(), endpt)[2]
                    
                    g2entry = grid.get_matrix_value(node.get_node_id(), endpt)[3]
                    
                    #print aentry
                    
                    
            
                    #h4[i] += c* lentry + g1entry * g1  + g2entry * g2 # G1, G2
                    h4[i] += c* lentry - g1entry * g1  - g2entry * g2 # -G1, -G2
                    
                    #print h4
                    
    return h4

h4 = h4(grid, Crhs, g1rhs, g2rhs)


def Lmatrix(grid):
    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
    
    for node in Nl:
        node.set_value(Nl.index(node))
#    
    
    
    Lmatrix = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt1)[0] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        Lmatrix[i, j] = lentry
    return Lmatrix
    

Lmatrix = Lmatrix(grid).todense()              

                
#Generate G1 matrix
def G1(grid):

    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
    
    for node in Nl:
        node.set_value(Nl.index(node))
    
    
    
    
    G1 = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
            #print node.get_coord(), i
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    #print  j
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    g1 = grid.get_matrix_value(node.get_node_id(), endpt1)[2] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        G1[i, j] = g1
    return G1

G1 = G1(grid).todense()

#def FD_G1(grid):
#    
#    Nl=[]
#    
#    for node in (not_slave_node(grid)):
#        
#        Nl.append([node,node.get_node_id()])
#        
#    Nl=sorted(Nl, key = itemgetter(1))
#    
#    Nl=[node[0] for node in Nl]
#    
#    h = 1/float(math.sqrt(len(Nl))+1)
#    
#    G1 = csr_matrix((len(Nl), len(Nl)))
#    
#    for node in Nl:
#        
#        node.set_value(Nl.index(node))
#        
#    for node in node_iterator(grid):
#        
#        if not node.get_slave():
#            
#            i = int(node.get_value())
#            
#            for endpt1 in connect_iterator(grid, node.get_node_id()):
#                
#                j = int(grid.get_value(endpt1))
#                
#                if not grid.get_slave(endpt1):
#                    
#                    if i-j == 0:
#                        
#                        G1[i,j] = -2*h
#                        
#                    elif abs(i -j) == np.sqrt(len(Nl)):
#                        
#                        #print i,j , '-1'
#                        
#                        G1[i,j] = h
#                        
#                        #print G1[i,j]
#        
#        
#    return G1
#
##
#G1 = FD_G1(grid).todense()

        
    
    




# Generate G2 matrix
    
def G2(grid):
    
    
    #build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)
    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
    
    for node in Nl:
        
        node.set_value(Nl.index(node))

    G2 = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    g2 = grid.get_matrix_value(node.get_node_id(), endpt1)[3] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        G2[i, j] = g2
                    
                    
    return G2



G2 = G2(grid).todense()

#def FD_G2(grid):
#    
#    Nl=[]
#    
#    for node in (not_slave_node(grid)):
#        
#        Nl.append([node,node.get_node_id()])
#        
#    Nl=sorted(Nl, key = itemgetter(1))
#    
#    Nl=[node[0] for node in Nl]
#    
#    h = 1/float(math.sqrt(len(Nl))+1)
#    
#    G2 = csr_matrix((len(Nl), len(Nl)))
#    
#    
#    
#    for node in Nl:
#        
#        node.set_value(Nl.index(node))
#        
#    for node in node_iterator(grid):
#        
#        if not node.get_slave():
#            
#            i = int(node.get_value())
#            
#            for endpt1 in connect_iterator(grid, node.get_node_id()):
#                
#                j = int(grid.get_value(endpt1))
#                
#                if not grid.get_slave(endpt1):
#                    
#                    if i-j == 0:
#                        
#                        G2[i,j] = -2*h
#                        
#                    elif abs(i -j) == 1:
#                        
#                        #print i,j , '-1'
#                        
#                        G2[i,j] = h
#                        
#                        #print G2[i,j]
#        
#        
#    return G2

#
#G2 = FD_G2(grid).todense()



                
                
# Generate dvector   

def dvector(grid, data):
    
    
    #build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  X, Y)
    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
    
    for node in Nl:
        node.set_value(Nl.index(node))
  


               
    dvector = zeros([len(Nl),1])
    
    for tri in Interior_triangle_iterator(grid):
        
        if tri[1].get_node_id() < tri[2].get_node_id():
            
            basis1 = set_polynomial_linear_2D(tri[0], tri[2], tri[1])
            
            Idi = int(tri[0].get_value())
            
            
            for i in range(len(Coord)):
                
                
                if NIn_triangle(tri[0] , tri[1], tri[2], Coord[i]):
                    
                        
                        
                        dvector[Idi,0] += basis1.eval(Coord[i][0], Coord[i][1]) * data[i]
                        
                        
                    
    return dvector/float(len(Coord))
#
dvector = dvector(grid, data)


rhs_d = dvector - h1
##
#
#
def Lets_Make_the_Damn_Big_Matrix():
    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
    
    ZeroMatrix = csr_matrix((len(Nl),len(Nl)))
    
    BigMat = bmat([[Amatrix, ZeroMatrix, ZeroMatrix, math.sqrt(alpha)*Lmatrix],\
                       [ZeroMatrix, Lmatrix, ZeroMatrix, -G1.T],\
                       [ZeroMatrix, ZeroMatrix, Lmatrix,  -G2.T],\
                       [math.sqrt(alpha)*Lmatrix, -G1, -G2, ZeroMatrix]])
    return BigMat
    
FEBigMat = Lets_Make_the_Damn_Big_Matrix().todense()
#print BigMat

Identity = np.identity(4*len(Nl))




#diag =np.diag(PD_BigMat)
#
#inv_diag = np.asarray([1/float(k) for k in diag])
#
#
#Diagonal = np.diag(inv_diag)
#
#Diagonal = np.sqrt(Diagonal)
#
#PrePD_BigMat = Diagonal* PD_BigMat * Diagonal







def plot_eigen(matrix, ind, nodenum, x1, y1):
    """ eigenvector corresponding to a given index, from big to small
    """
    
    eigenvalue, eigenvector = LA.eig(matrix)
    print 'max', max(Leigenvalue), 'min', min(Leigenvalue), 'cond', LA.cond(Lmatrix)
    
    eigvalue = abs(Leigenvalue)
    eigvalue = eigvalue.tolist()
    
    # Find the index of the smallest eigvalue
    Aeigvalue = eigvalue.index(sorted(eigvalue, rever = True)[ind])
    
    Aeigvector = eigenvector[:, Aeigvalue]
    Aeigvector = Aeigvector[0: nodenum]
    Aeigvector = np.reshape(Aeigvector, (np.sqrt(nodenum), np.sqrt(nodenum)))
    
    ax = Axes3D(plt.gcf())
    ax.plot_surface(x1,y1, Aeigvector)
    
    
        #w eigenvalues, v eigenvectors
    Leigenvalue,Leigvector=LA.eigh(FEPD_BigMat)
    #eigvector, eigvalue, vh=np.linalg.svd(Amatrix)
    print 'max', max(Leigenvalue), 'min', min(Leigenvalue), 'cond', LA.cond(FEPD_BigMat)
    
    
    Leigvalue=abs(Leigenvalue)
    Leigvalue=Leigvalue.tolist()
    
    
    #print  'min' , min(abs(eigvalue))
    
    # Return the indices of the max value along an axis.
    
    Lmaxeigvalue = Leigvalue.index(sorted(Leigvalue, reverse= True)[-1])
    
    # Find the eigenvector corresponding to max eigenvalue
    Lmaxeigvector = Leigvector[:,Lmaxeigvalue]
    
    
    # Take the real component 
    
    
    #Take the section of eigenvector for c
    Lmaxeigvector=Lmaxeigvector[len(Nl):2*len(Nl)]
    
    
    Lmaxeigvector = np.reshape(Lmaxeigvector, (int(np.sqrt(len(Nl))), int(np.sqrt(len(Nl)))))
    
    
    
    
    ax = Axes3D(plt.gcf())
    
    ax.plot_surface(x1,y1, Lmaxeigvector)
    
    pass
    



rhs = zeros([4*len(Nl),1])
rhs[0:len(Nl)]=rhs_d
rhs[len(Nl):2*len(Nl)]=-h2/float(math.sqrt(alpha))
rhs[2*len(Nl):3*len(Nl)]=-h3/float(math.sqrt(alpha))
rhs[3*len(Nl):4*len(Nl)]=-math.sqrt(alpha)*h4

value_vector =spsolve(FEBigMat, rhs)
error_vector = copy(value_vector)[0:len(Nl)]





for node in node_iterator(grid):
    
    # The error is only calculated at the interior nodes
    
    if not node.get_slave():
        
        # What row in the matrix does the current node correspond to
        i = int(node.get_value())
        

        value = value_vector[i]
        coord = node.get_coord()
        error_vector[i] = error_vector[i]-l3(coord)
        #print value
        
        #print value

        
        # Record the value at the current node
        node.set_value(value)
        
        
        
    # If the node is a slave node, its value should be given by the boundary condition
    else:
        node.set_value(Cube(node.get_coord()[0],node.get_coord()[1]))
        

print "Error norm :", linalg.norm(error_vector)*h
    
  

def plot_fem_solution(grid):
    from grid.Grid import Grid
    from grid.NodeTable import node_iterator
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from matplotlib.pyplot import show, figure
    from scipy.interpolate import griddata
    from numpy import mgrid, array
#   

   
# The connection values are set separately with respect to the location, 

    #grid=Grid()
#    
    #Find the position of the nodes and the values
    node_x=[]
    node_y=[]
    node_v=[]
    for node in node_iterator(grid):
        coord = node.get_coord()
        node_x.append(coord[0])
        node_y.append(coord[1])
        node_v.append(node.get_value())
        
    # Store the results in an array
    
    node_x = array(node_x)
    node_y = array(node_y)
    node_v = array(node_v)
    
    #print('node_x',node_x)
    #print('node_value', node_v)
    
    
    # Initialise the figure
    fig = plt.figure(figsize=(8,5)) 
    ax = fig.gca(projection='3d') 
    
    
    # Interpolate the nodes onto a structured mesh
#    X, Y = mgrid[node_x.min():node_x.max():10j,
#                 node_y.min():node_y.max():10j]
    
    X, Y = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    
    Z = griddata((node_x,node_y), node_v, (X,Y), method='cubic')
    
    
    # Make a surface plot
    ax.plot_surface(X, Y, Z,cmap='viridis',
                       linewidth=0)
#    surf=ax.plot_surface(X, Y, Z,cmap='viridis',
#                       linewidth=0)
    
    # Set the z axis limits
    #ax.set_zlim(node_v.min(),node_v.max())
    
    # Make the ticks looks pretty
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    
    # Include a colour bar
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    # Show the plot
    show()

plot_fem_solution(grid)











