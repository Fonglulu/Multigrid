#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 15:44:05 2018

@author: FENG Shi Lu
"""
import math
from copy import copy
from grid.GlobalNodeID import GlobalNodeID
from grid.Node import Node
from grid.function.FunctionStore import zero, exp_soln
from grid.Edge import Edge, DomainSet, RefineSet
    
    
#######################################################################
# add_edge
#Add a new edge on table. Give the attributes to added edge joined by two points, 
#and for the added edge we store from both directions
#######################################################################

def add_edge(grid, id1, id2, refine_type, location, \
        boundary_function = exp_soln):

    edge1 = Edge(id1, id2) #Edge Class
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
# add_Packman_edges
# Add the edge in two nodes on grid. For boundary and interior nodes,
# we have separated cases. 
#
#######################################################################
    
    
    
def add_Packman_edges(grid, N_length, N_angle, node_mesh,boundary_function):
    base = RefineSet.base_edge
    not_base = RefineSet.not_base_edge
    intp = DomainSet.interior
    bnd = DomainSet.boundary
    
    #create the 0 node   
    coord= [0,0]
    glboal_id=GlobalNodeID()
    #glboal_id.set_no(mesh_integers_2D(i,j))
    glboal_id.set_level(0)
    node=grid.create_node(glboal_id,coord,True,0.0)
    #print(i-1,j)
    grid.add_node(node)
    zeroid = node.get_node_id()
    
       #nodes at the 0 
    
    id1 = node_mesh[0][0].get_node_id()
    add_edge(grid, zeroid, id1, not_base, bnd, boundary_function)
    
    id1  = node_mesh[0][-1].get_node_id()
    add_edge(grid, zeroid, id1, not_base, bnd, boundary_function)
    
    for i in range(1, N_angle-1):
        id1 =node_mesh[0][i].get_node_id()
        add_edge(grid, zeroid, id1, not_base, intp)
    
    
    
    for i in range(0, N_length-2):
        for j in range(1, N_angle-1):
            #print(i,j)
            id1 = node_mesh[i][j].get_node_id()
            id2 = node_mesh[i+1][2*j-1].get_node_id()
            add_edge(grid, id1, id2, not_base, intp)
            
            id2 = node_mesh[i+1][2*j+1].get_node_id()
            add_edge(grid, id1, id2, not_base, intp)
            
            id2 = node_mesh[i+1][2*j].get_node_id()
            add_edge(grid, id1, id2, not_base, intp)
            
            id2 = node_mesh[i][j+1].get_node_id()
            add_edge(grid, id1, id2, base, intp)
            
            
        N_angle=2*N_angle-1    
        print(N_angle)
        
    #nodes on  upper stright boundary, N_length=5
    for i in range(0, N_length-2):
        id1 = node_mesh[i][0].get_node_id()
        id2 = node_mesh[i+1][0].get_node_id()
        add_edge(grid, id1, id2, not_base, bnd, boundary_function)
         
        id2 = node_mesh[i+1][1].get_node_id()
        add_edge(grid, id1, id2, not_base, intp)
        
        id2 = node_mesh[i][1].get_node_id()
        add_edge(grid, id1, id2, base, intp) 
        
    #nodes on lower stright boundary 
    for i in range(0, N_length-2):
        id1 = node_mesh[i][-1].get_node_id()
        id2 = node_mesh[i+1][-1].get_node_id()
        add_edge(grid, id1, id2, not_base, bnd, boundary_function)
        
        id2 = node_mesh[i+1][-2].get_node_id()
        add_edge(grid, id1, id2, not_base,intp)
        

    #nodes on curve boundary
    for i in range(1, N_angle):
        id1 = node_mesh[-1][i].get_node_id()
        
        id2 = node_mesh[-1][i-1].get_node_id()
        add_edge(grid, id1, id2, base, bnd, boundary_function)
        
 

        
        


#######################################################################
# add_square_nodes
# By knowing the radius and angle, use polar coordinate to give coordinate
# of noedes. Assign the coordinates and NodeID to the node_mesh created by
# generate_node_mesh
#######################################################################
def add_Packman_nodes(grid, N_length, N_angle, node_mesh, start_angle, end_angle):
    from grid.GlobalNodeID import mesh_integers_2D
    
    angle=end_angle-start_angle
    
    for i in range(1, N_length):
        for j in range(N_angle):
            radius = (i+1) / float(N_length)
            theta=(angle /( N_angle-1)) * j
            
            coord= [radius * math.cos(theta), radius * math.sin(theta)]
            
            glboal_id=GlobalNodeID()
            glboal_id.set_no(mesh_integers_2D(i,j))
            glboal_id.set_level(0)
            node=grid.create_node(glboal_id,coord,False,0.0)
        
            node_mesh[i-1][j] = copy(node)
            grid.add_node(node)

        N_angle = 2*N_angle-1
        

    
    
    
    for i in range(N_length-2):
        #print(i)
        node=node_mesh[i][0]
        grid.set_slave(node.get_node_id(),True)
        node=node_mesh[i][-1]
        grid.set_slave(node.get_node_id(),True)
        
    for i in range(len(node_mesh[-1])):
        node=node_mesh[-1][i]
        grid.set_slave(node.get_node_id(),True)
        
    
 

        

#######################################################################
# generate_node_mesh
# This function generates the node_mesh with the strcuture such that the 
# first level has 5 nodes, the second has 9 nodes so on and so forth.
#######################################################################
        
        
def generate_node_mesh(N_angle, N_length, grid):
    
    node=Node()
    node_mesh=[]
    for i in range(1, N_length):
        node_level=[]
        for j in range(N_angle):
            node_level.append(node)
        N_angle=2*N_angle - 1
        node_mesh.append(node_level)
    return node_mesh







#######################################################################
# build_packman_grid
# Given certian amount of levels and nodes on first level, along with 
# the angle, this function creates the nodes and builds the edges.
#######################################################################

def build_packman_grid(boundary_function, start_angle , end_angle , N_angle, N_length, grid):

        
    node_mesh = generate_node_mesh(N_angle,N_length,grid)

    add_Packman_nodes(grid, N_length, N_angle, node_mesh, start_angle, end_angle) 
    
    add_Packman_edges(grid, N_length, N_angle, node_mesh,boundary_function)    
        

