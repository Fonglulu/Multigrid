#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 09:10:28 2018

@author: shilu
"""

import numpy as np
from grid.NodeTable import node_iterator
from grid.ConnectTable import connect_iterator 
from grid.function.FunctionStore import exp_soln, exp_rhs, sin_soln, sin_rhs
from DDinit import DD_init
from copy import copy
from scipy import eye, zeros, linalg
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import splu, spsolve


    
    
    
    
    
    
    
    
    
    
    
    
############################################################################    
# DD_preconditioner
#
# p = D^{-1}*r, where D is a block diagonal matrix and block i corresponds 
# to the part of the fem matrix stored with worker i. DD_init should be called 
# before calling this routine as it finds the sparese LU decomposition of 
# the FEM natrix. Nodes should not be added or removed from the grid in between
# calling DD_init and DD_precondtioner
#
#
# The entries in ghost_vector_table are not used. It is listed in the set of
# parameters fro consistency.
#
# It is assumed that r is the residual vector so the values at the slave of 
# node is zero.
#
#
#Input: Grid grid
#     : dictionary ghost_vector_table
#     : dictionary vector_table
#     : int r_index
#     : int p_index
#
# Output: THe result is stroed in the p_index of the arrays
# stored in vector_table. The slvae node entries are not changed
############################################################################
def DD_preconditioner(grid, ghost_vector_table, vector_table, r_index, 
                      p_index):
    
    r=[]
    i = 0
    for node in node_iterator(grid):
        node_id = node.get_node_id()
        #print (node_id._id_no,'number')
        if not node.get_slave():
            
            #print(node.get_coord())
            
            vector = vector_table[str(node_id)]
            
            #print(node_id._id_no, vector[r_index], 'vv')
            
            r.append(vector[r_index])
            #print('first loop:',i)
            i+=1

    #matrix = spsolve(DD_init(grid),r)
    matrix = DD_init(grid).solve(np.array(r))
    #vector[p_index] = vector[r_index]
    #vector_table[str(node_id)] = vector
    j=0
    for node in node_iterator(grid):
        node_id = node.get_node_id()
        
        if not node.get_slave():
            vector = vector_table[str(node_id)]
            vector[p_index] = matrix[j]
            vector_table[str(node_id)] = vector
            #print('second loop:',j)
            j+=1
            
            
            
            
            
            
            
            
            
            
            
            
#            vector = vector_table[str(node_id)]
#                    
#            vector[p_index] = vector[r_index]
#                    
#            vector_table[str(node_id)] = vector
#            
#            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
    
