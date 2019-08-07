#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 15:29:04 2018

@author: shilu
"""

from grid.Grid import Grid
from grid.function.FunctionStore import zero
from BuildSquare import build_square_grid

# Comment out the following when imported into ghost grid

grid= Grid()
n=5
build_square_grid(n,grid,zero)




 #######################################################################
# plot_fem_grid
# print the list of triangle and  plot the grid
#######################################################################
def plot_fem_grid(c, grid):
    from grid.NodeTable import node_iterator 
    from grid.NghNodes import ngh_iterator
    from grid.EdgeTable import endpt_iterator 
    from grid.Edge import DomainSet
    from matplotlib.pyplot import plot, show
    import matplotlib.pyplot as plt
    
    #for ngh in grid.reference_ghost_table():
        #print(ngh.get_coord(),'ghost')
    for node in node_iterator(grid):
        print(node.get_coord(),'full nodes')
        #print (node.get_global_id()._id_no, node.get_coord())
            #print(node)
            #plot(node._coord)
        for endpt in endpt_iterator(grid,node.get_node_id()):
                #print(endpt, grid.get_coord(endpt))
                if grid.is_in(endpt):
                    #print(grid.get_coord(endpt), 'm')
                    x = [node.get_coord()[0], grid.get_coord(endpt)[0]]
                    y = [node.get_coord()[1], grid.get_coord(endpt)[1]]
                    #plot(x,y)
                if  grid.reference_ghost_table().is_in(endpt):
                          print(grid.reference_ghost_table().get_coord(endpt),'ghost')
                          x = [node.get_coord()[0], grid.reference_ghost_table().get_coord(endpt)[0]]
                             #print(ngh[1][i]._id_no)
                             #print(grid.get_coord(ngh[1][i]._id_no))
                          y = [node.get_coord()[1], grid.reference_ghost_table().get_coord(endpt)[1]]
                 
       
                plot(x,y)   
                #print(endpt._id_no,"connected")
                #plot(node._coord,grid.get_node(endpt._id_no).get_coord())
    for ghost in ngh_iterator(grid.reference_ghost_commun()):
        #print(ghost[1])
        for node in ghost[1]:
            #print(grid.reference_ghost_table().get_coord(node),'coord')
            for endpt in endpt_iterator(grid, node):
                #print(endpt)
                if  grid.reference_ghost_table().is_in(endpt):
                    x = [grid.reference_ghost_table().get_coord(node)[0], grid.reference_ghost_table().get_coord(endpt)[0]]
                    y = [grid.reference_ghost_table().get_coord(node)[1], grid.reference_ghost_table().get_coord(endpt)[1]]
                    plot(x,y)
            #for endpt in endpt_iterator()
    plt.title('subgrid'+str(c))          
    show ()
#    
#    
#plot_fem_grid(grid)

        
    # print the list of triangles
#    for node in node_iterator(grid):
#        if not node.get_slave():
#            #print(node)
#        #node_i=node.get_node_id()
#        node_i=[node.get_node_id()._id_no]
#        #print(node_i)
#        for endpt_1 in endpt_iterator(grid,node.get_node_id()):
#            node_j=[endpt_1._id_no]
#            #print(node_j)
#            for endpt_2 in endpt_iterator(grid, node.get_node_id()):
#                node_k=[endpt_2._id_no]
                #print(node_k)
                #if grid.is_edge(endpt_1,endpt_2) == True:
                    #print(node_i+node_j+node_k)
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
