#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 13:21:28 2018

@author: Feng
"""

from grid.Grid import Grid
from grid.function.FunctionStore import zero
from BuildPackman import build_packman_grid
from BuildSquare import build_square_grid
import matplotlib.pyplot as plt
import numpy as np

#grid= Grid()
#N_length=4
#N_angle=4
#start_angle=0
##The input for end_angle has to be float
#end_angle=1.0
#
#
#build_packman_grid(zero, start_angle, end_angle, N_angle, N_length, grid)

#grid= Grid()
#n=4
#build_square_grid(n,grid,zero)


#######################################################################
# plot_fem_grid:
# Give triangle and plot the grid
#######################################################################
def plot_fem_grid(grid):
    from grid.NodeTable import node_iterator 
    from grid.EdgeTable import endpt_iterator 
    from grid.Edge import DomainSet
    from matplotlib.pyplot import plot, show
    
 #This prints out the list of triangles
    for node in node_iterator(grid):
        #node_i=node.get_node_id()
        node_i=[node.get_node_id()._id_no]
        #print(node_i)
        for endpt_1 in endpt_iterator(grid,node.get_node_id()):
            node_j=[endpt_1._id_no]
            #print(node_j)
            for endpt_2 in endpt_iterator(grid, node.get_node_id()):
                node_k=[endpt_2._id_no]
                #print(node_k)
                #if grid.is_edge(endpt_1,endpt_2) == True:
                   # print(node_i+node_j+node_k)
    
  

# This plot the grid out from build_packman_grid    
    for node in node_iterator(grid):
        
        for endpt in endpt_iterator(grid, node.get_node_id()):
            x = [node.get_coord()[0], grid.get_node(endpt._id_no).get_coord()[0]]
            y = [node.get_coord()[1], grid.get_node(endpt._id_no).get_coord()[1]]
  
            plot(x,y)
    x1 = np.linspace(0+1/float(12), 1.0-1/float(12),12)
    y1 = np.linspace(0+1/float(12), 1.0-1.0/float(10),12)
    X, Y = np.meshgrid(x1,y1)
    #fig, ax = plt.subplots()
    #fig,ax = plt.subplots(1,1,figsize=(10,10))
    plt.scatter(X, Y,s =1)
    show ()
    

    

