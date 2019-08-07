# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This module stores the routine that subdivides the grid into the given
number of sections. It firstly uses metis to find a non-overlapping 
division of the full nodes. Given the full nodes, edges and connections, it
termines which nodes should be added as ghost nodes for each subdivision.
Finally it builds the, full and ghost, neighbour node table for each
subdivision. 

Note, the routines assumes the edge and connection tables have been 
defined.
"""

# import appropriate information from other classes

from copy import deepcopy
from Grid import Grid
from pymetis import part_graph
from NodeTable import NodeTable, node_iterator
from EdgeTable import endpt_iterator
from ConnectTable import connect_iterator
    
#######################################################################
# add_connections
#
# Loop through the full nodes in sub_grid and copy across the connection
# stars from grid.
#
# Input: Grid sub_grid
#        Grid grid
#
# Output: The all of the connection stars corresponding to the full nodes
# in sub_grid have been copied across from grid. The matrix values are
# also copied across. 
#
#######################################################################

def add_connections(sub_grid, grid):
    """Copy the connection stars from grid into sub_grid"""
    
    # Loop over the full nodes in sub_grid
    for node in node_iterator(sub_grid):
        node_id = node.get_node_id()
        
        # Loop over the connection star stored in grid
        for endpt in connect_iterator(grid, node_id):
            
            # Copy the connection across
            stiff = grid.get_matrix_value(node_id, endpt)
            sub_grid.add_connection(node_id, endpt, stiff)
                
            
    
    

#######################################################################
# add_edges
#
# Loop through the full nodes in sub_grid and copy across the edge
# stars from grid.
#
# Input: Grid sub_grid
#        Grid grid
#
# Output: The all of the edge stars corresponding to the full nodes
# in sub_grid have been copied across from grid.
#
#######################################################################

def add_edges(sub_grid, grid):
    """Copy the edge stars from grid into sub_grid"""
    
    # Loop over the full nodes in sub_grid
    for node in node_iterator(sub_grid):
        node_id = node.get_node_id()
        
        # Loop over the connection star stored in grid
        for endpt in endpt_iterator(grid, node_id):

            # Copy the edge across
            edge = grid.get_edge(node_id, endpt)
            sub_grid.add_edge_all(edge)
            
    
    
#######################################################################
# add_ghost_nodes
#
# The routine loops over the edge and connections stars in sub_grid and
# if any of the endpts are not in the full node table they are added to
# the ghost node table
#
#
# Input: Grid sub_grid
#        Grid grid
#
# Output: if (node_id1, node_id2) is an edge or a connection and 
# node_id2 is not in the full node table, it is added as a ghost node
#
#######################################################################
def add_ghost_nodes(sub_grid, grid):
    """Add the ghost nodes to complete the connections to the full nodes"""
    
    # Find the ghost node table
    ghost_table = sub_grid.reference_ghost_table()
    
    # Loop over the full nodes in the grid
    for node in node_iterator(sub_grid):
        node_id = node.get_node_id()
        
        # Loop over the edge stars
        for endpt in endpt_iterator(grid, node_id):
           
            # If the endpoint is not alread in sub_grid, add it as a ghost node
            if not sub_grid.is_in(endpt) and not ghost_table.is_in(endpt):
                node = grid.get_node(endpt)
                ghost_table.add_node(node)
                
        # Loop over the connection stars
        for endpt in connect_iterator(grid, node_id):
            
            # If the endpoint is not alread in sub_grid, add it as a ghost node
            if not sub_grid.is_in(endpt) and not ghost_table.is_in(endpt):
                node = grid.get_node(endpt)
                ghost_table.add_node(node)
            

#######################################################################
# add_ghost_node_stars
#
# The routine loops the ghost nodes in the sub_grid and adds an edge 
# (or connection) to another full or ghost node in the sub_grid if there
# is a corresponing edge (or connection) the original grid
#
# The routine assumes that the ghost node table has been added to the
# subgrid
#
# Input: Grid sub_grid
#        Grid grid
#
# Output: Suppose node_id1 is a ghost node in sub_grid and node_id2 is
# either a full node or a ghost, then if (node_id1, node_id2) is an 
# edge in grid it is added to sub_grid. Similarly if (node_id1, node_id2) 
# is a connection to grid it is added to sub_grid, and the matrix value
# is copied across.
#
#######################################################################
def add_ghost_node_stars(sub_grid, grid):
    """Add the edges/connections from the ghost nodes"""
    
    # Find the ghost node table
    ghost_table = sub_grid.reference_ghost_table()
    
    # Loop over the ghost nodes in sub_grid
    for node in node_iterator(ghost_table):
        node_id = node.get_node_id()
        
        # Loop over the edge star in grid
        for endpt in endpt_iterator(grid, node_id):
            
            # If the endpoint is in sub_grid, add the edge to sub_grid
            if sub_grid.is_in(endpt) or ghost_table.is_in(endpt):
               sub_grid.add_edge_all(grid.get_edge(node_id, endpt))
                               
        # Loop over the connection stars in grid
        for endpt in connect_iterator(grid, node_id):
            
            # If the endpoint is in sub_grid, add the connection to sub_grid
            if sub_grid.is_in(endpt) or ghost_table.is_in(endpt):
                sub_grid.add_connection(node_id, endpt)
                stiff = grid.get_matrix_value(node_id, endpt)
                sub_grid.set_matrix_value(node_id, endpt, stiff)
                

                
                
#######################################################################
# add_full_nodes
#
# Copy the full nodes from grid into sub_grid. full_node_partition is
# the node partioning given by metis.
#
#
# Input: integer worker_no 
#        Grid sub_grid
#        Grid grid
#        Adjacency list full_node_partition
#
# Output: The edge between id1 and id2 as well as the edge 
# from id2 to id1 has been added to grid
#
#######################################################################
def add_full_nodes(worker_no, sub_grid, grid, full_node_partition):
    """Copy the full nodes from grid into sub_grid"""
    
    # Get the list of nodes to be added to the current worker
    full_node_list = full_node_partition[worker_no]
    
    # Copy the nodes across from grid into sub_grid
    for node_id in full_node_list:
        node = grid.get_node(node_id)
        sub_grid.add_node(node)
            
#######################################################################
# build_sub_grid
#
# Build the sub_grid corresponding to worker_no. This includes the 
# full nodes, ghost nodes, edges and connections. It does not include
# the neighbour node table, that is done by build_commun.
#
#
# Input: integer worker_no 
#        Grid sub_grid
#        Grid grid
#        Adjacency list full_node_partition
#
# Output: The full node table, ghost node table, edge table and 
# connection table for the subgrid to be stored on worker_no is returned
# in the grid.
#
#######################################################################
def build_sub_grid(worker_no, grid, full_node_partition):
    """Build a subgrid"""
    
    # Initialise the grid
    sub_grid = Grid()
    
    # Add the full nodes to the grid
    add_full_nodes(worker_no, sub_grid, grid, full_node_partition)
    
    # Add the ghost nodes, by completing the edges and connections
    add_ghost_nodes(sub_grid, grid)
    
    # Copy the edge stars across (for the full nodes)
    add_edges(sub_grid, grid)
    
    # Copy the connection stars across (for the ghost nodes)
    add_connections(sub_grid, grid)
    
    # Add in the necessary edges and connections going from the ghost nodes
    add_ghost_node_stars(sub_grid, grid)
    
    # Return the subgrid
    return sub_grid
  

#######################################################################
# build_commun
#
# Build the neighbour node tables for all of the subgrids. 
#
# It is assumed that the full and ghost node table for each subgrid
# has been built.
#
#
# Input: dictionary sub_grids
#        Grid grid
#        Adjacency list full_node_partition
#
# Output: The full and ghost neighbour node table for each sub_grid
# has been built
#
#######################################################################
def build_commun(sub_grids, grid, full_node_partition):
    """Build the neighbour node tables"""
    
    # We use this temporary table to find which worker will have the 
    # full node copy
    tmp_node_table = NodeTable()
    for worker in full_node_partition:
        for node_id in full_node_partition[worker]:
            node = grid.get_node(node_id)
            node.set_value(worker)
            tmp_node_table.add_node(node)
    
    # Loop over the subgrids
    for worker in sub_grids:
        
        # Work on the current subgrid
        sub_grid = sub_grids[worker][0]

        # Find the ghost table corresponding to the current subgrid
        ghost_table = sub_grid.reference_ghost_table()
        
        # Find the ghost neighbour node table corresponding to the current 
        # subgrid
        ghost_commun = sub_grid.reference_ghost_commun()

        # Loop over all of the nodes in the ghost table
        for node in ghost_table:
            node_id = node.get_node_id()
            
            # Find which worker contains the full node copy
            full_worker = tmp_node_table.get_value(node_id)
            
            # Let the worker with the full node copy know that this worker
            # contains a ghost node copy
            full_commun = sub_grids[full_worker][0].reference_full_commun()
            full_commun.add_neighbour_node(worker, node_id)

            # Keep a record of which worker contains the full node copy of the
            # current ghost node
            ghost_commun.add_neighbour_node(full_worker, node_id)
           
        
        
#######################################################################
# build_sub_grids
#
# Given the grid paritioning in full_node_partition, build all of the
# subgrids
#
#
# Input: Grid grid
#        Adjacency list full_node_partition
#
# Output: A dictionary the same size as full_node_partition is created
# The key is the worker number and the value is the subgrid that should
# be sent to that worker_no
#######################################################################
def build_sub_grids(grid, full_node_partition):
    """Build all of the subgrids"""
    
    # create an empty dictionary
    sub_grids = {}
    
    # Loop over all of the workers
    for worker in full_node_partition:
        
        # Build the subgrid for the given worker no
        sub_grid = build_sub_grid(worker, grid, full_node_partition)
        
        # Store the subgrid in the dictionary
        sub_grids.setdefault(worker, []).append(deepcopy(sub_grid))
        
    
    # Given the set of sub_grids, find the communication pattern
    build_commun(sub_grids, grid, full_node_partition)
    
    return sub_grids
    
#######################################################################
# divide_full_nodes
#
# Use metis to divide the full nodes up into a non-overlapping set of
# sub partions. It used the edges to form an adjaceny list.
#
#
# Input: integer no_workers
#        Grid grid
#
# Output: The nodes and edges defining the square domain as in the grid
#
#######################################################################
def divide_full_nodes(no_workers, grid):
    """Partition the grid into the given number of workers"""

    # Initialise a list of nodes
    node_list = list()
    
    # Initialise the dictionary to hold the node partitions
    full_partition = {}
    
    # Initialise the adjacency list used by metis to find the node partitions
    adjacency = {}
    
    # Build a list of nodes (this is primarily done to assign each node an
    # integer value)
    for node in node_iterator(grid):
        node_list.append(node.get_node_id())
            
    # Build the adjacency list
    
    # Loop over the nodes
    for node in node_iterator(grid):
        node_id = node.get_node_id()
        node_index = node_list.index(node_id)
        
        # Loop over the edges
        for endpt in endpt_iterator(grid, node_id):
            endpt_index = node_list.index(endpt)
            
            # Make a note of the relation between the two endpoints
            adjacency.setdefault(node_index, []).append(endpt_index)
            
    # Use metis to partition the grid
    c, partition = part_graph(no_workers, adjacency)
    
    # Make note of which node should go to which worker
    for i in range(len(partition)):
        full_partition.setdefault(partition[i], []).append(node_list[i])
        
    # Return the partitioning
    return full_partition
 
#######################################################################
# subdivie_grid
#
# Subdivide the grid into no_workers subgrids.
#
#
# Input: integer no_workers
#        Grid grid
#
# Output: A dictionary the same size as no_workers is created.
# The key is the worker number and the value is the subgrid that should
# be sent to that worker_no
#
####################################################################### 
def subdivide_grid(no_workers, grid):
    
    full_node_partition = divide_full_nodes(no_workers, grid)
    
    sub_grids = build_sub_grids(grid, full_node_partition)
    return sub_grids