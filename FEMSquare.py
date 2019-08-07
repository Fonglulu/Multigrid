#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 09:00:06 2018

@author: shilu
"""

# Import appropriate information from other classes
from grid.Grid import Grid
from BuildSquare import build_square_grid
#from FemPlot import plot_fem_solution
from grid.NodeTable import node_iterator
from grid.ConnectTable import connect_iterator   
from grid.function.FunctionStore import exp_soln, exp_rhs, sin_soln, sin_rhs
from copy import copy
from scipy import eye, zeros, linalg
import numpy as np
from numpy import inf

from BuildEquation import build_equation_linear_2D, Poisson_tri_integrate, build_matrix_fem_2D

    
# Initialise the grid
grid = Grid()

# Build a 5*5 grid
i = 2
n = 2**i

x = np.linspace(0,1.1,11)
y = np.linspace(0,1.1,11)

X, Y = np.meshgrid(x,y)
# Specify the boundary conditions and the rhs
true_soln = sin_soln
rhs = sin_rhs

# Create a square FEM grid of size n*n
build_square_grid(n, grid, true_soln)

#build_equation_linear_2D(grid, Poisson_tri_integrate, rhs)
build_matrix_fem_2D(grid, Poisson_tri_integrate, X, Y, rhs)


# Loop through the grid and find all of the node in the interior of the
# domain. At the same time make a record of which row in the matrix 
# corresponds to which node
count = 0
for node in node_iterator(grid):
    if not node.get_slave():
        node.set_value(count)
        count = count+1

# Initialise the A matrix
fem_matrix = eye(count)

# Initialise the load vector
load_vector = zeros([count, 1])

# Define the A matrix and load vector
for node in node_iterator(grid):
    
    # Ignore slave (or boundary) nodes
    if not node.get_slave():
        
        # Which row corresponds to the current node?
        i = int(node.get_value())
        
        # Add in the entry to the load vector
        coord = node.get_coord()
        load_vector[i] = grid.get_load(node.get_node_id())#load_vector[i]+rhs(coord)
        
        # Loop over the matrix entries in the current row
        for endpt1 in connect_iterator(grid, node.get_node_id()):
            
            # Which column corresponds to the current node?
            j = int(grid.get_value(endpt1))
            
            # What is the corresponding matrix value (in the FEM grid)
            stiffness = grid.get_matrix_value(node.get_node_id(), endpt1)[1]
            
            #print stiffness
            
            # We must not include slave nodes in the matrix columns
            if not grid.get_slave(endpt1):
                fem_matrix[i, j] = stiffness
                
                #print fem_matrix[i,j]
                
            # Update the load vector to take non-zero boundary conditions
            # into account
            else:
                coord = grid.get_coord(endpt1)
                load_vector[i] =load_vector[i]-stiffness*true_soln(coord)
 
# Solve the system Ax = b
value_vector = linalg.solve(fem_matrix, load_vector)

# Initialise a vector to keep track of the errors

error_vector = copy(value_vector)


# Store the solution back in the FEM grid and calcuate the error at
# each node
for node in node_iterator(grid):
    
    # The error is only calculated at the interior nodes
    
    if not node.get_slave():
        
        # What row in the matrix does the current node correspond to
        i = int(node.get_value())
        
        # Find the error at the current node
        coord = node.get_coord()
        value = value_vector[i]
        error_vector[i] = error_vector[i]-true_soln(coord)
        
        # Record the value at the current node
        node.set_value(value)
        
    # If the node is a slave node, its value should be given by the
    # boundary condition
    else:
        node.set_value(true_soln(node.get_coord()))
        


# Calculate the norm of the error
print "Error norm :", linalg.norm(error_vector, inf)

# Plot the solution
#plot_fem_solution(grid)
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