# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This files contains routines to create a finite element grid defined on
a square domain. The main routines are build_square_grid and 
build_square_grid_matrix.

build_square_grid builds an n*n uniform finite element grid defined on 
a square domain (0, 1)* (0, 1) with the given Dirichlet boundary conditions.

build_square_grid_matrix is similar to build_square_grid but it also uses
the five points stencil to construct the stiffness matrix that is defined
on the square domain. The stiffness matrix is stored in the connection table.


"""

# Import appropriate information from other classes
from copy import copy
from grid.GlobalNodeID import GlobalNodeID
from grid.Node import Node
from grid.function.FunctionStore import zero, linear,sqr_rhs
from grid.Edge import Edge, DomainSet, RefineSet
import numpy as np
    
#######################################################################
# add_edge
#
# Add an edge between id1 and id2 to the FEM grid. The edge location,
# boundary function and refinement type is also set.
#
#
# Input: Grid grid
#        NodeID id1 
#        NodeID id2
#        RefineSet refine_type
#        DomainSet location
#        Function boundar_function
#
# Output: The edge between id1 and id2 as well as the edge 
# from id2 to id1 has been added to grid
#
#######################################################################

def add_edge(grid, id1, id2, refine_type, location, \
        boundary_function=sqr_rhs):
    """ Add an edge to the grid""" 
    edge1 = Edge(id1, id2)
    edge1.set_location(location)
    edge1.set_refine_type(refine_type)
    edge1.set_boundary_function(boundary_function)
    edge1.set_refine_level(0)
    edge2 = copy(edge1)
    edge2.set_endpt1(id2)
    edge2.set_endpt2(id1)
    grid.add_edge_all(edge1)
    grid.add_edge_all(edge2)
    
    
#######################################################################
# add_square_connections
#
# Given the n*n mesh of nodes in a square grid, add all of the connections 
# needed to store the 5 point stencil matrix. The connections 
# corresponding to the diagonal entries in the grid are given a value of
# 0. This is not necessary, but is what will happen when the matrix is
# calculated automatically.
#
# Input: Grid grid
#        Integer n
#        Array of Nodes node_mesh
#
# Output: The connection between id1 and id2 as well as the connection 
# from id2 to id1 has been added to grid. The values assigned are 
# given by the 5-point matrix stencil.
#
#######################################################################

def add_square_connections(grid, n, node_mesh):
    """ Add all of the connections in a square grid""" 
    
    h = 1.0/(n-1)
    h2 = h*h
    
    # Loop through all of the edges in the interior of the domain
    for i in range(0, n-1):
        for j in range(0, n-1):
            
            # Add the diagonal matrix entry
            id1 = node_mesh[i][j].get_node_id()
            grid.add_connection(id1, id1, [4.0/h2])
            
            # Add a connection from the current node to the one above it
            id2 = node_mesh[i+1][j].get_node_id()
            grid.add_connection(id1, id2, [-1.0/h2])
            grid.add_connection(id2, id1, [-1.0/h2])
           
            # Add a connection from the current node to the one to the right
            id2 = node_mesh[i][j+1].get_node_id()
            grid.add_connection(id1, id2, [-1.0/h2])
            grid.add_connection(id2, id1, [-1.0/h2])
            
            # Add a connection from the current node to the one diagonally above
            id2 = node_mesh[i+1][j+1].get_node_id()
            grid.add_connection(id1, id2, [0.0])
            grid.add_connection(id2, id1, [0.0])

            
    # Now work on the edges around the boundary of the domain
            
    # Loop over the edges along the right boundary
    for i in range(0, n-1):
        
        # Add the diagonal matrix entry
        id1 = node_mesh[n-1][i].get_node_id()
        grid.add_connection(id1, id1, [4.0/h2])

        
        # Add an edge from the current node to the one above it
        id2 = node_mesh[n-1][i+1].get_node_id()
        grid.add_connection(id1, id2, [-1.0/h2])
        grid.add_connection(id2, id1, [-1.0/h2])
         
            
    # Loop over the edges along the top boundary
    for i in range(0, n-1):
        
        # Add the diagonal matrix entry
        id1 = node_mesh[i][n-1].get_node_id()
        grid.add_connection(id1, id1, [4.0/h2])

        
        # Add an edge from the current node to the one to the right
        id2 = node_mesh[i+1][n-1].get_node_id()
        grid.add_connection(id1, id2, [-1.0/h2])
        grid.add_connection(id2, id1, [-1.0/h2])
 
            
    
    
#######################################################################
# add_square_edges
#
# Given the n*n mesh of nodes in a square grid, add all of the edges needed
# to form a collection of triangles. The routine assumes that nodes in the
# mesh have been assigned their local id.
#
#
# Input: Grid grid
#        Integer n
#        Array of Nodes node_mesh
#        Function boundar_function
#
# Output: The edge between id1 and id2 as well as the edge 
# from id2 to id1 has been added to grid
#
#######################################################################

def add_square_edges(grid, n, node_mesh, boundary_function):
    """ Add all of the edges in a square grid""" 
    
    # For short hand notation
    base = RefineSet.base_edge
    not_base = RefineSet.not_base_edge
    intp = DomainSet.interior
    bnd = DomainSet.boundary
    
    # Loop through all of the edges in the interior of the domain
    for i in range(1, n-1):
        for j in range(1, n-1):
            
            # Add an edge from the current node to the one above it
            id1 = node_mesh[i][j].get_node_id()
            id2 = node_mesh[i+1][j].get_node_id()
            add_edge(grid, id1, id2, not_base, intp)
            
            # Add an edge from the current node to the one to the right
            id2 = node_mesh[i][j+1].get_node_id()
            add_edge(grid, id1, id2, not_base, intp)
            
            # Add an edge from the current node to the one diagonally above
            id2 = node_mesh[i+1][j+1].get_node_id()
            add_edge(grid, id1, id2, base, intp)

            
    # Now work on the edges around the boundary of the domain
    
    # Loop over the edges along the left boundary (ignore the nodes at the end)
    for i in range(1, n-1):
        
        # Add an edge from the current node to the one to the right
        id1 = node_mesh[0][i].get_node_id()
        id2 = node_mesh[1][i].get_node_id()
        add_edge(grid, id1, id2, not_base, intp, boundary_function)
        
        # Add an edge from the current node to the one above it
        id2 = node_mesh[0][i+1].get_node_id()
        add_edge(grid, id1, id2, not_base, bnd, boundary_function)
        
        # Add an edge from the current node to the one diagonally above
        id2 = node_mesh[1][i+1].get_node_id()
        add_edge(grid, id1, id2, base, intp, boundary_function)
            
    # Loop over the edges along the right boundary
    for i in range(0, n-1):
        
        # Add an edge from the current node to the one above it
        id1 = node_mesh[n-1][i].get_node_id()
        id2 = node_mesh[n-1][i+1].get_node_id()
        add_edge(grid, id1, id2, not_base, bnd, boundary_function)
            
    # Loop over the edges along the bottom boundary 
    # (ignore the nodes at the end)
    for i in range(1, n-1):
        
        # Add an edge from the current node to the one above it  
        id1 = node_mesh[i][0].get_node_id()
        id2 = node_mesh[i][1].get_node_id()
        add_edge(grid, id1, id2, not_base, intp, boundary_function)
        
        # Add an edge from the current node to the one to the right
        id2 = node_mesh[i+1][0].get_node_id()
        add_edge(grid, id1, id2, not_base, bnd, boundary_function)
        
        # Add an edge from the current node to the one diagonally above
        id2 = node_mesh[i+1][1].get_node_id()
        add_edge(grid, id1, id2, base, intp, boundary_function)
            
    # Loop over the edges along the top boundary
    for i in range(0, n-1):
        
        # Add an edge from the current node to the one to the right
        id1 = node_mesh[i][n-1].get_node_id()
        id2 = node_mesh[i+1][n-1].get_node_id()
        add_edge(grid, id1, id2, not_base, bnd, boundary_function)
            
    # The previous boundary loops miss the edges going from the bottom left
    # corner
    
    # Add an edge from the corner to the one to above it
    id1 = node_mesh[0][0].get_node_id()
    id2 = node_mesh[1][0].get_node_id()
    add_edge(grid, id1, id2, not_base, bnd, boundary_function)
    
    # Add an edge from the corner to the one to the right
    id2 = node_mesh[0][1].get_node_id()
    add_edge(grid, id1, id2, not_base, bnd, boundary_function)
    
    # Add an edge from the corner to the one diagonally above
    id2 = node_mesh[1][1].get_node_id()
    add_edge(grid, id1, id2, base, intp, boundary_function)    


#######################################################################
# add_square_nodes
#
# Given the mesh of nodes in a square grid, add them to the finite 
# element grid. The means that both the local and gloval ids are 
# assigned to the nodes. The routine also assigns a local id to the 
# nodes in the mesh.
#
# It is assumed that the domain is (0,1)x(0,1)
#
#
# Input: Grid grid
#        Integer n
#        Array of Nodes node_mesh
#        Function boundary_function
#
# Output: The edge between id1 and id2 as well as the edge 
# from id2 to id1 has been added to grid
#
#######################################################################
def add_square_nodes(grid, n, node_mesh):
    """ Add all of the nodes in a square grid""" 
    
    # Import appropriate information from other modules
    from numpy import arange
    from grid.GlobalNodeID import mesh_integers_2D

    # Grid space
    h = 1.0/(n-1)
    spacing = arange(0.0, 1.0+h, h)
    
    # Loop through the nodes in the mesh and find their coords as well as
    # local and global ids
    
    for i in range(n):
        x = spacing[i]
        for j in range(n):
            y = spacing[j]
            
            # Find the node coordinate
            coord = [x, y]
            
            # Find the global id
            global_id = GlobalNodeID()
            global_id.set_no(mesh_integers_2D(i, j))
            global_id.set_level(0)
            
            # Create a new node
            node = grid.create_node(global_id, coord, False, 0.0)
            
            # Add the node to the grid
            grid.add_node(node)
            
            # Record the lodal node id in the node mesh
            node_mesh[i][j] = copy(node)
            
    # Set the boundary node to be slave nodes
    for i in range(n):
        node = node_mesh[0][i]
        grid.set_slave(node.get_node_id(), True)
        node= node_mesh[n-1][i]
        grid.set_slave(node.get_node_id(), True)
        node = node_mesh[i][0]
        grid.set_slave(node.get_node_id(), True)
        node = node_mesh[i][n-1]
        grid.set_slave(node.get_node_id(), True)
        

#######################################################################
# build_square_grid
#
# Build a square finite element grid of size n*n, with the given 
# Dirichlet boundary conditions. n should be > 2.
#
# It is assumed that the domain is (0,1)x(0,1)
#
# Input: integer n
#        Grid grid
#        Function Dirichelt
#
# Output: The nodes and edges defining the square domain as in the grid
#
#######################################################################
def build_square_grid(n, grid, boundary_function):
    """Build a square finite element grid"""
    
    # Initialsie a node
    node = Node()
    
    # This is a mesh of nodes
    node_mesh = [[node]*n for x in xrange(n)]
    
    # Assume there are at least 3 nodes in each direction
    assert n > 2, "the grid size must be > 2"
            
    # Add the nodes to the fem grid and find their local and global ids
    add_square_nodes(grid, n, node_mesh)
             
    # Add the edges that join the nodes together 
    add_square_edges(grid, n, node_mesh, boundary_function)
        

#######################################################################
# build_square_grid_matrix
#
# Build a square finite element grid of size n*n, with the given 
# Dirichlet boundary conditions. The 5 point stencil matrix is also
# stored with the grid. n should be > 2.
#
# It is assumed that the domain is (0,1)x(0,1)
#
# Input: integer n
#        Grid grid
#        Function Dirichelt
#
# Output: The nodes and edges defining the square domain as in the grid
#
#######################################################################
def build_square_grid_matrix(n, grid, boundary_function):
    """Build a square finite element grid as well as the matrix"""
    
    # Initialsie a node
    node = Node()
    
    # This is a mesh of nodes
    node_mesh = [[node]*n for x in xrange(n)]
    
    # Assume there are at least 3 nodes in each direction
    assert n > 2, "the grid size must be > 2"
            
    # Add the nodes to the fem grid and find their local and global ids
    add_square_nodes(grid, n, node_mesh)
             
    # Add the edges that join the nodes together 
    add_square_edges(grid, n, node_mesh, boundary_function)    

    # Add square connections
    add_square_connections(grid, n, node_mesh)
    
