# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This class defines an iterator that extracts the triangles in a grid
"""

# Import appropriate information from other classes
from grid.NodeTable import node_iterator
from grid.EdgeTable import endpt_iterator

    
#######################################################################
# triangle_iterator
#
# Use the node_iterator and endpt_iterator generators to define an
# iterator that loops over the triangles in the grid. 
#
# Note that the iterator does not check the order. That is it will return
# (node1, node2, node3), (node1, node3, node2), .....
#
#
# Input: Grid grid
#
# Output: a list containing the three nodes sitting on the vertices of
# the triangle
#
#######################################################################
def triangle_iterator(grid):
    """Iterate over the triangles in a grid"""
    
    # Get the ghost node table
    ghost_table = grid.reference_ghost_table()

    # Loop over the nodes in the grid
    for node in node_iterator(grid):
        
        # Loop over the edges joined to the node
        for endpt1 in endpt_iterator(grid, node.get_node_id()): 
            
            # Get the node sitting to the endpoint (which may be a
            # ghost node)
            if grid.is_in(endpt1):
                node1 = grid.get_node(endpt1)
            else:
                node1 = ghost_table.get_node(endpt1)
               
            
            #Loop over another set of edges joined to the node
            for endpt2 in endpt_iterator(grid, node.get_node_id()):
                
                # If the endpoints are joined by an additional edge
                if grid.is_edge(endpt1, endpt2):
                    
                    # Get the node sitting to the endpoint (which may be a
                    # ghost node)
                    if grid.is_in(endpt2):
                        node2 = grid.get_node(endpt2)
                    else:
                        node2 = ghost_table.get_node(endpt2)
 
                    # Then return the triangle
                    yield [node, node1, node2]
        
                      
                    
                    
def Interior_triangle_iterator(grid):
    
    """Iterate over the interior triangles in a grid"""
    
    # Get the ghost node table
#    ghost_table = grid.reference_ghost_table()

    # Loop over the nodes in the grid
    for node in node_iterator(grid):
        
        if not node.get_slave():
            
            # Loop over the edges joined to the node
            for endpt1 in endpt_iterator(grid, node.get_node_id()): 
                
                # Get the node sitting to the endpoint (which may be a
                # ghost node)
                
                if grid.is_in(endpt1):
                    node1 = grid.get_node(endpt1)
#                else:
#                    node1 = ghost_table.get_node(endpt1)
                   
                
                #Loop over another set of edges joined to the node
                for endpt2 in endpt_iterator(grid, node.get_node_id()):
                    
                    # If the endpoints are joined by an additional edge
                    if grid.is_edge(endpt1, endpt2):
                    
                        # Get the node sitting to the endpoint (which may be a
                        # ghost node)
                        if grid.is_in(endpt2):
                            node2 = grid.get_node(endpt2)
#                        else:
#                            node2 = ghost_table.get_node(endpt2)
     
                        # Then return the triangle
                        yield [node, node1, node2]
        
    
