#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 13:17:37 2018

@author: shilu
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-


##############################################################################
    ##########   ##      ##   ###########     ###########       ###########
        ##       ##      ##       ##          ##                         ##
        ##       ##      ##       ##          ##                         ##
        ##       ##      ##       ##          ##                         ##
        ##       ##      ##       ##          ##########        ########### 
        ##       ##      ##       ##          ##                ##
        ##       ##      ##       ##          ##                ##
        ##        ########        ##          ###########       ########### 
# This works
##############################################################################


# import appropriate information from other classes
from copy import copy
import numpy as np
from grid.Grid import Grid
from BuildSquare import build_square_grid
from BuildPackman import  build_packman_grid
from BuildEquation import build_equation_linear_2D, Poisson_tri_integrate
from grid.NodeTable import node_iterator
from grid.EdgeTable import  endpt_iterator
from grid.ConnectTable import connect_iterator
from grid.function.FunctionStore import exp_soln, exp_rhs, sin_soln, sin_rhs,zero
from scipy import eye, zeros, linalg
from numpy import inf
from operator import itemgetter
from scipy.sparse import csc_matrix,csr_matrix
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


# Initialise the grid
grid = Grid()


# Build a 5*5 grid
i = 3
n= 2**i+1

# Specify the boundary contiditions and the rhs

true_soln =sin_soln
rhs = sin_rhs

# Create a square FEM grid od size n*n
build_square_grid(n, grid, rhs)

#Initialise the A matrix
fem_matrix = zeros([81,81])

# Initialise the load vector
load_vector = zeros([49,1])


build_equation_linear_2D(grid, Poisson_tri_integrate, rhs)


def not_slave_node(grid):
    for node in node_iterator(grid):
        if not node.get_slave():
            yield node

####################################################################
            ####### Creat Load Vector ######
IntNode=[]
for node in (not_slave_node(grid)):
    #print(i)
    IntNode.append([node,node.get_node_id()])
IntNode_ID=sorted(IntNode, key = itemgetter(1))
IntNode_Ordered=[node[0] for node in IntNode_ID]
#print(IntNode_Ordered,IntNode_Ordered[0].get_node_id()._id_no)


for i in range(len(IntNode_Ordered)):
    load_value=IntNode_Ordered[i].get_load()
    #print(load_value)
    load_vector[i]= load_vector[i]+load_value
    
    
for j, node in enumerate(IntNode_Ordered):
        boundary_sum=[]
        id1=node.get_node_id()
        coord=node.get_coord()
        for endpt in (connect_iterator(grid, id1)):
            if grid.get_slave(endpt):
                connect_value=grid.get_matrix_value(id1, endpt)
                #print(node.get_node_id()._id_no,endpt._id_no)
                boundary_sum.append(true_soln(grid.get_coord(endpt))*connect_value)
        load_vector[j]=load_vector[j]-sum(boundary_sum)
        #print(load_vector)
        
        
        
####################################################################       
        
        
# Enumerate interior NodeID in order of ID
NodeID_Load=[]
for q, node in enumerate(not_slave_node(grid)):
    NodeID_Load.append(node.get_node_id())
NodeID_Load.sort()
#print(NodeID_Load)



 #For every NodeID in NodeID_load, find its connecting nodes, reorder them and assign them in stiffness matrix
#for i, nodeid in enumerate(NodeID_Load):
#    Endpt=[]
#    for j, endpt1 in enumerate(connect_iterator(grid, nodeid)):
#        if not grid.get_slave(endpt1):
#            Endpt.append(endpt1)
#    Endpt.sort() 
#    print(Endpt)
#    #B=[]
#    for k,endpt2 in enumerate(Endpt):
#        #print(endpt2)
#        #B.append([k,endpt2])
#        stiffness = grid.get_matrix_value(nodeid,endpt2)
#        print(i,k,nodeid._id_no,endpt2._id_no,stiffness)
#        fem_matrix[nodeid._id_no-6,endpt2._id_no-6]=stiffness
##def not_slave_node(grid):
##    for node in node_iterator(grid):
##        if not node.get_slave():
#            yield node
#
#

# For interior nodes and their interior connections, assign them in stiffness matrix (25*25) accordingly
for node in (NodeID_Load):
    Endpt=[]
    for endpt in (connect_iterator(grid, node)):
        if not grid.get_slave(endpt):
            Endpt.append(endpt)
        Endpt.sort()
    for endpt2 in (Endpt):
        stiffness = grid.get_matrix_value(node,endpt2)
        #print(node._id_no,endpt2._id_no,stiffness)
        fem_matrix[node._id_no,endpt2._id_no]=stiffness

matrixfem =np.matrix(fem_matrix)
reduced_matrix=matrixfem[~(matrixfem==0).all(1).A1]
transfem=np.transpose(reduced_matrix) 
redumatrix=transfem[~(transfem==0).all(1).A1]           
        
             
             
BoundaryNode = []
for i, node in enumerate(node_iterator(grid)):
    if node.get_slave():
        BoundaryNode.append([node,node.get_node_id()])
BoundaryNodeLoad_ID = sorted(BoundaryNode,key = itemgetter(1))
BoundaryNodeLoad_Ordered=[node[0] for node in BoundaryNodeLoad_ID]

             
boundary_value=[]
for node in (BoundaryNodeLoad_Ordered):
    boundary_value.append(true_soln(node.get_coord()))



#print(load_vector)
#   
#
#
###### L =[ (q, NODE_ID(q))], items list from small to large  #####
#L=[]
#for q, node in enumerate(node_iterator(grid)):
#    #print(q)
#    L.append(node.get_node_id())
#L.sort()
#
#
#
####### i is the nodes that has been sorted in order, for each node i, sort the 
####### the connected endpoint in order, put it in list B
#for i, nodeid in enumerate(L):
#    A=[]
#    for j, endpt1 in enumerate(connect_iterator(grid, nodeid)):
#        A.append(endpt1)
#    A.sort() 
#    #B=[]
#    for k,endpt2 in enumerate(A):
#        #B.append([k,endpt2])
#        stiffness = grid.get_matrix_value(nodeid,endpt2)
#        fem_matrix[nodeid._id_no,endpt2._id_no]=stiffness
#        print(i,k,nodeid._id_no,endpt2._id_no,stiffness)
#        #print(nodeid._id_no,endpt2._id_no,stiffness)
#        
##print(fem_matrix)
#
value_vector = linalg.solve(redumatrix , load_vector)
       
        
#print(fem_matrix)               
# Solve the system Ax = b
#value_vector = linalg.solve(fem_matrix, load_vector)


true_vector=[]
for i in range(len(IntNode_Ordered)):
    true_vector.append(true_soln(IntNode_Ordered[i].get_coord()))
    np.array(true_vector)   
#print(fem_matrix)               
# Solve the system Ax = b
#value_vector = linalg.solve(fem_matrix, load_vector)


 #Initalise a vector to keep track of the errors
error_vector = copy(value_vector)-copy(true_vector)




# Calculate the norm of the error
print "Error nor:",  linalg.norm(error_vector, inf)



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
    for i in range(len(NodeID_Load)):
        grid.set_value(NodeID_Load[i]._id_no, value_vector[i])
    
    for i in range(len(BoundaryNodeLoad_Ordered)):
        grid.set_value(BoundaryNodeLoad_Ordered[i].get_node_id()._id_no, boundary_value[i])
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
    
    print('node_x',node_x)
    print('node_value', node_v)
    
    
    # Initialise the figure
    fig = figure()
    ax = fig.gca(projection='3d') 
    ax = fig.gca()
    
    
    # Interpolate the nodes onto a structured mesh
    X, Y = mgrid[node_x.min():node_x.max():10j,
                 node_y.min():node_y.max():10j]
    
    Z = griddata((node_x,node_y), node_v, (X,Y), method='cubic')
    
    
    # Make a surface plot
    surf = ax.plot_surface(X, Y, Z,
            cmap=cm.coolwarm , linewidth=0, antialiased=False)
    
    # Set the z axis limits
    ax.set_zlim(node_v.min(),node_v.max())
    
    # Make the ticks looks pretty
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    
    # Include a colour bar
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    # Show the plot
    show()

        
        
plot_fem_solution(grid)     
#        
#        
        
        