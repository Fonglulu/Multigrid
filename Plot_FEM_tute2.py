#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 15:40:44 2018

@author: FENG Shi Lu
"""

##############################################################################
    ##########   ##      ##   ###########     ###########       ###########
        ##       ##      ##       ##          ##                         ##
        ##       ##      ##       ##          ##                         ##
        ##       ##      ##       ##          ##                         ##
        ##       ##      ##       ##          ##########        ########### 
        ##       ##      ##       ##          ##                ##
        ##       ##      ##       ##          ##                ##
        ##        ########        ##          ###########       ########### 
##############################################################################


# import appropriate information from other classes
from copy import copy
import numpy as np
from grid.Grid import Grid
from BuildSquare import build_square_grid_matrix
from grid.NodeTable import node_iterator
from grid.ConnectTable import connect_iterator
from grid.EdgeTable import  endpt_iterator
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
i = 2
n= 2**i+1

# Specify the boundary contiditions and the rhs
#true_soln = zeros
#rhs = sin_rhs

true_soln =sin_soln
rhs = sin_rhs

# Create a square FEM grid od size n*n
build_square_grid_matrix(n, grid, true_soln)


#Initialise the A matrix
fem_matrix = eye(49)

#fem_matrix = lil_matrix((225,225))

# Initialise the load vector
load_vector = zeros([49,1])


def not_slave_node(grid):
    for node in node_iterator(grid):
        if not node.get_slave():
            yield node



#ALoad = []
#for i, node in enumerate(node_iterator(grid)):
#    ALoad.append([node,node.get_node_id()])
#ALoad_node_ID = sorted(ALoad,key = itemgetter(1))
#ALoad_node=[node[0] for node in ALoad_node_ID]
#    


############################################################
# Generate Boundary nodes and list them in order of ID and find the boundary values
            
BoundaryNode = []
for i, node in enumerate(node_iterator(grid)):
    if node.get_slave():
        BoundaryNode.append([node,node.get_node_id()])
BoundaryNodeLoad_ID = sorted(BoundaryNode,key = itemgetter(1))
BoundaryNodeLoad_Ordered=[node[0] for node in BoundaryNodeLoad_ID]
    
    
boundary_value=[]
for node in (BoundaryNodeLoad_Ordered):
    boundary_value.append(true_soln(node.get_coord()))
    
############################################################





############################################################
#  Generates Interior nodes and list them in order of ID
    
IntNode=[]
for i, node in enumerate(not_slave_node(grid)):
    #print(i)
    IntNode.append([node,node.get_node_id()])
IntNode_ID=sorted(IntNode, key = itemgetter(1))
IntNode_Ordered=[node[0] for node in IntNode_ID]
############################################################



############################################################

 #Add boundary value from corresponding load vector, if it is zero boundary condition, it will add nothing
for j, node in enumerate(IntNode_Ordered):
        boundary_sum=[]
        id1=node.get_node_id()
        coord=node.get_coord()
        for endpt in endpt_iterator(grid, id1):
            if grid.get_slave(endpt):
                if abs(node.get_node_id()._id_no-endpt._id_no)==1 or \
                   abs(node.get_node_id()._id_no-endpt._id_no)==9:
                       #print(node.get_node_id()._id_no,endpt._id_no)
                       boundary_sum.append(true_soln(grid.get_coord(endpt)))
        load_vector[j]=load_vector[j]+sum(boundary_sum)*64


        
# Find the load vector unfer zero Dirichlet boundary condition        
for i in range(len(IntNode_Ordered)):
    coord=IntNode_Ordered[i].get_coord()
    load_vector[i]= load_vector[i]+rhs(coord)
###############################################################
    
#print(load_vector)
   

##############################################################
# Enumerate interior NodeID in order of ID
NodeID_Load=[]
for q, node in enumerate(not_slave_node(grid)):
    #print(q)
    NodeID_Load.append(node.get_node_id())
NodeID_Load.sort()



# For every NodeID in NodeID_load, find its connecting nodes, reorder them and assign them in stiffness matrix
for i, nodeid in enumerate(NodeID_Load):
    Endpt=[]
    for j, endpt1 in enumerate(connect_iterator(grid, nodeid)):
        Endpt.append(endpt1)
    Endpt.sort() 
    #B=[]
    for k,endpt2 in enumerate(Endpt):
        #B.append([k,endpt2])
        stiffness = grid.get_matrix_value(nodeid,endpt2)
        fem_matrix[i,k]=stiffness*float(64)
        #print(i,k,nodeid._id_no,endpt2._id_no,stiffness)
        #print(nodeid._id_no,endpt2._id_no,stiffness)

##############################################################        
#print(fem_matrix)


#############################################################
        #Solve for the apprxoimating nodal value and find the true value at those nodes
value_vector = linalg.solve(fem_matrix , load_vector)
       
true_vector=[]
for i in range(len(IntNode_Ordered)):
    true_vector.append(true_soln(IntNode_Ordered[i].get_coord()))
    np.array(true_vector)
    
#print(fem_matrix)               
# Solve the system Ax = b
#value_vector = linalg.solve(fem_matrix, load_vector)


# Initalise a vector to keep track of the errors
error_vector = copy(value_vector)-copy(true_vector)
#
#
#
#
## Calculate the norm of the error
print "Error norm:", linalg.norm(error_vector, inf)


#        
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
        
        
        
        
        
        
        
        
        
        
        
        
        