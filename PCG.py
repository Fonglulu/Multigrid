# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This file contains the preconditioned conjugate gradient routine. 

"""

# Import appropriate information from other classes
from grid.NodeTable import node_iterator
from grid.ConnectTable import connect_iterator
from math import sqrt
    
#######################################################################
# matrix_mult
#
# Evaluate q = Ad. The A matrix is stored in the grid. The values for
# d are stored in ghost_vector_table and vector_table and are indexed
# by d_index. The values for q are returned in vector_table and are 
# index by q_index. Note that grid and ghost_vector_table are 
# read but are not updated. 
#
# As this is just intended to be a helper routine, no checking is done.
# That means, for example, the full nodes stored in grid have an 
# equivalent entry in vector_table and the ghost nodes stored in grid 
# have an equivalent entry in ghost_vector_table.
#
# Input: Grid grid
#      : dictionary ghost_vector_table
#      : dictionary vector_table
#      : int d_index
#      : int q_index
#
# Output: The result is stored in the q_index index of the arrays
# stored in vector_table. The slave node entries are not changed.
#
#######################################################################     
def matrix_mult(grid, ghost_vector_table, vector_table, d_index, q_index):
    """ Evaluate q = Ad"""

    # Loop over the full nodes
    for node in node_iterator(grid):
        node_id = node.get_node_id()
        
        # Do not change the slave nodes
        if not node.get_slave():
            
            # Find sum A_{ij} d_j
            sum = 0.0
            
            # Loop over the end points connect to node i
            for endpt in connect_iterator(grid, node_id):
                
                # Find (A)_{ij}
                stiffness = grid.get_matrix_value(node_id, endpt)
                
                # Find d_j, which may be a full or ghost node
                if grid.is_in(endpt):
                    d = vector_table[str(endpt)][d_index]
                else:
                    d = ghost_vector_table[str(endpt)][d_index]
                    
                # Update the sum
                sum = sum + stiffness*d
                
            # Store the result
            vector_table[str(node_id)][q_index] = sum


#######################################################################
# calculate_residual
#
# Evaluate r = b-Ax. The A matrix is stored in the grid. The values for
# x are stored in ghost_vector_table and vector_table and are indexed
# by x_index. The values for b and r are returned in vector_table and are 
# index by b_index and r_index. Note that grid and ghost_vector_table are 
# read but are not updated. 
#
# As this is just intended to be a helper routine, no checking is done.
# That means, for example, the full nodes stored in grid have an 
# equivalent entry in vector_table and the ghost nodes stored in grid 
# have an equivalent entry in ghost_vector_table.
#
# Input: Grid grid
#      : dictionary ghost_vector_table
#      : dictionary vector_table
#      : int x_index
#      : int b_index
#      : int r_index
#
# Output: The result is stored in the r_index index of the arrays
# stored in vector_table. The slave node entries are not changed.
#
#######################################################################     
def calculate_residual(grid, ghost_vector_table, vector_table,
                       x_index, b_index, r_index):
    """ Evaluate r = b-Ax"""

    
    # Find Ax
    matrix_mult(grid, ghost_vector_table, vector_table, 
                x_index, r_index)
                
    # Find r = b - Ax
    
    # Loop over the full nodes
    for node in node_iterator(grid):
        node_id = node.get_node_id()
        
        # Ignore slave nodes
        if not node.get_slave():
            
            # Find r_i = b_i - (Ax)_i
            vector = vector_table[str(node_id)]
            vector[r_index] = vector[b_index]-vector[r_index]
            vector_table[str(node_id)] = vector
            

#######################################################################
# inner_produce
#
# Evaluate <v1, v2>. The vector v1 is indexed by index1 in vector_table,
# the vector v2 is indexed by index2 in vector_table.
#
# As this is just intended to be a helper routine, no checking is done.
# That means, for example, index1 and index2 are well defined indices
#
# Input: dictionary vector_table
#      : int index1
#      : int index2
#
# Output: a float number is returned. The sum is taken over all of the
# entries in vector_table
#
#######################################################################          
def inner_product(vector_table, index1, index2):
    """Evaluate <v1, v2>"""
    
    # Initialise the sum
    sum = 0.0
    
    # Find sum v1_i*v2_i
    for key, value  in vector_table.iteritems():
        sum = sum+value[index1]*value[index2]
        
    # Return sum v1_i*v2_i
    return sum

#######################################################################
# conjugate_gradient
#
# Solve the system Ax = b using the preconditioned conjugate gradient
# method. A is stored in the connection table. x is stored in the value
# field of the nodes in the grid. b is stored in the load field of the 
# nodes in the grid.
#
# Input: Communication commun
#      : Grid grid
#      : function preconditioner
#      : real tolerance
#      : int max_no_loops
#
# Output: The result is stored in the r_index index of the arrays
# stored in vector_table. The slave node entries are not changed.
#
#######################################################################     
def conjugate_gradient(commun, grid, preconditioner, 
                       tolerance = 1.0E-12, max_no_loops = 5000):
                           
    """Preconditioned conjugate gradient method"""
    
    # These are indices into temporary storage areas needed by the routine
    residual_index = 0
    direction_index = 1
    tmp_index = 2
    precond_index = 3
    x_index = 4
    b_index = 5
    
    # Initialise the scaling terms used by CG
    resid_norm2 = 0.0
    alpha_new = 0.0
    alpha_old = 0.0
    beta = 0.0
    gamma = 0.0
    
    # Explicitly recalculate the residual after this number of steps
    residual_recalc = 50
    
    # How many loops did the CG method take to solve the problem?
    no_loops = max_no_loops
    
    
    # With each full node in the grid, store an array containing all of the
    # information needed by the CG method
    
    # We use a dictionary to store the information
    vector_table = dict()
    
    # Loop of the full nodes
    for node in node_iterator(grid):
        
        # 6 pieces of information must be stored with each node
        vector = [0.0]*6
        
        # Find the initial guess
        vector[x_index] = node.get_value()
        
        # Find the right hand side
        vector[b_index] = node.get_load()
        
        # Store the array
        vector_table[str(node.get_node_id())] = vector
        
    # With each ghost node in the grid, store an array containing all of the
    # information needed by the CG method    
    ghost_table = grid.reference_ghost_table()
    
    # We use a dictionary to store the information
    ghost_vector_table = dict()
    
    # Loop over the ghost nodes
    for node in node_iterator(ghost_table):
        
        # 6 pieces of information must be stored with each node
        vector = [0.0]*6
        
        # Store the array
        ghost_vector_table[str(node.get_node_id())] = vector
        
    # Update the ghost nodes
    commun.update_vector_table(grid, vector_table, ghost_vector_table)
    
    # Find the initial residual
    calculate_residual(grid, ghost_vector_table, vector_table,
                       x_index, b_index, residual_index)
    
    # Residual norm^2 = <r, r>
    resid_norm2 = inner_product(vector_table, residual_index, residual_index)
    resid_norm2 = commun.global_sum(resid_norm2)
    
    # Apply the preconditioner (s = M^{-1}r)
    preconditioner(grid, ghost_vector_table, vector_table, residual_index, 
                   precond_index)
                       
    # Initialise the direction vector (d = s)
    for key, value in vector_table.iteritems():
        value[direction_index] = value[precond_index]
       
    # Find alpha_new = < r, s>
    alpha_new = inner_product(vector_table, residual_index, precond_index)
    alpha_new = commun.global_sum(alpha_new)

    # Start the CG iterations
    for loop in range(max_no_loops):
        
        # Temporarily store A*direction_vector (q = A*d)
        commun.update_vector_table(grid, vector_table, ghost_vector_table)
        matrix_mult(grid, ghost_vector_table, vector_table, direction_index, 
                    tmp_index)
                    
        # Beta = alpha_new/<d, q>
        dq = inner_product(vector_table, direction_index, tmp_index)
        dq = commun.global_sum(dq)
        beta = alpha_new/dq
                   
        # Update the current solution in the given direction (x = x + beta*d)
        for key, value  in vector_table.iteritems():
            value[x_index] = value[x_index] + beta*value[direction_index]
        
        # Find the current residual
        
        # Every residual_recalc iterations recalculate the residual explicitly
        # (r = b-A*x)
        if loop%residual_recalc == 0:
            commun.update_vector_table(grid, vector_table, ghost_vector_table)
            calculate_residual(grid, ghost_vector_table, vector_table, 
                               x_index, b_index, residual_index)
        else:
            
        # Otherwise use quick formula to update the residual (r = r-beta*q)
            for key, value in vector_table.iteritems():
                value[residual_index] = value[residual_index] \
                - beta*value[tmp_index]
 
        # Residual norm^2 = <r, r>
        resid_norm2 = inner_product(vector_table, residual_index, 
                                    residual_index)
        resid_norm2 = commun.global_sum(resid_norm2)

        # If the residual is small, then exit
        if (sqrt(resid_norm2) < tolerance):
            no_loops = loop
            break
        
        # alpha_old = alpha_new
        alpha_old = alpha_new
        
        # Apply the preconditioner (s = M^{-1}r)
        commun.update_vector_table(grid, vector_table, ghost_vector_table)
        preconditioner(grid, ghost_vector_table, vector_table, residual_index, 
                       precond_index)

        # Find alpha_new = < r, s>
        alpha_new = inner_product(vector_table, residual_index, precond_index)
        alpha_new = commun.global_sum(alpha_new)
        
        # Gamma = alpha_new/alpha_old
        gamma = alpha_new/alpha_old
        
        # Update the direction vector (d = r + gamma d)
        for key, value in vector_table.iteritems():
            value[direction_index] = value[precond_index] \
            + gamma*value[direction_index]

    # Transfer the results from the temporary data structure into the grid
    # data structure
    for node in node_iterator(grid):
        node.set_value(vector_table[str(node.get_node_id())][x_index])
        
    # Return the number of CG iterations as well as the residual
    return no_loops, sqrt(resid_norm2)