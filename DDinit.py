#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 19:10:10 2018

@author: shilu
"""

############################################################################
# DD_int
#
# Find the LU decompostion of the matrix stored in this processor
#
# DD_init must be called, once, before calling DD_preconditioner
#
#Input: Grid grid
#
#Output: The LU deomposition is found and stored in a global variable
#
############################################################################
import numpy as np
from grid.NodeTable import node_iterator
from grid.ConnectTable import connect_iterator 
from grid.function.FunctionStore import exp_soln, exp_rhs, sin_soln, sin_rhs
from copy import copy
from scipy import eye, zeros, linalg
from scipy.sparse import lil_matrix, csc_matrix, csr_matrix
from scipy.sparse.linalg import splu, spsolve, LinearOperator

def DD_init(grid):
    
    
        count = 0
        for node in node_iterator(grid):
             if not node.get_slave():
                  #print(node.get_coord())
                  node.set_value(count)
                  count = count+1
        #print count, 'count'
        #print(count,'count')
        # Initialise the A matrix
        #fem_matrix = csr_matrix((count, count), dtype=float)
        fem_matrix = eye(count)
        
        # Initialise the load vector
        #load_vector = zeros([count, 1])
        
        # Define the A matrix and load vector
        for node in node_iterator(grid):
            
            # Ignore slave (or boundary) nodes
            if not node.get_slave():
                
                # Which row corresponds to the current node?
                i = int(node.get_value())
                
                # Add in the entry to the load vector
                #coord = node.get_coord()
                #load_vector[i] = load_vector[i]+rhs(coord)
                
                # Loop over the matrix entries in the current row
                for endpt1 in connect_iterator(grid, node.get_node_id()):

                    if grid.is_in(endpt1):
                        
                        j = int(grid.get_value(endpt1))
                        
                        stiffness = grid.get_matrix_value(node.get_node_id(), endpt1)
                        
                        if not grid.get_slave(endpt1):  
                            
                            fem_matrix[i, j] = stiffness
                    
#                    if grid.reference_ghost_table().is_in(endpt1):
#                        
#                        j = int(grid.reference_ghost_table().get_value(endpt1))
#                        
#                        #stiffness = grid.get_matrix_value(node.get_node_id(), endpt1)
#                        
#                        if not grid.reference_ghost_table().get_slave(endpt1):  
#                            
#                            fem_matrix[i, j] = 0
                    
                    #if not grid.is_in(endpt1):
                        #print "not in grid"
                    # Which column corresponds to the current node?
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    #stiffness = grid.get_matrix_value(node.get_node_id(), endpt1)  
                    
                    # We must not include slave nodes in the matrix columns
                    #if not grid.get_slave(endpt1):
                        #fem_matrix[i, j] = stiffness
        #print(fem_matrix)
                        
                    # Update the load vector to take non-zero boundary conditions
                    # into account
                    #else:
                        #coord = grid.get_coord(endpt1)
                        #load_vector[i] = load_vector[i]-stiffness*true_soln(coord)
       
        #print fem_matrix
#        fem= splu(fem_matrix)
#        fem = LinearOperator(fem.shape,matvec=fem.solve)
        return splu(fem_matrix)
        #return fem_matrix
       
