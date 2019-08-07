#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 19:08:41 2018

@author: shilu
"""

###########################################################################
# identity_preconditioner
# p = I*r. That is, the contents index by r_index are copied into the part
# of the array indexed by p_index.




from grid.NodeTable import node_iterator




def identity_preconditioner(grid, ghost_vector_table, vector_table, 
                            r_index, p_index):
    """ p =I *r"""
    
    
    # Loop through the list of full nodes
    for node in node_iterator(grid):
        node_id = node.get_node_id()
        
        # Ignore the slace nodes
        if not node.get_slave():
            
            # Get the array of values assigned to the current node
            vector = vector_table[str(node_id)]
           
            
            # Copy across the entries indexed by r_index into those
            # indexed by p_index
            vector[p_index] = vector[r_index]
            
            # Update the array of values assignmed to the current node
            vector_table[str(node_id)] = vector
            
            
    #print('...............')

