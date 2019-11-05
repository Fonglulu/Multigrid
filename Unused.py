#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 08:20:25 2019

@author: shilu
"""

import numpy as np
#######################################################################################
# This routine is a collection of unused rountines
#######################################################################################
def h1_bd(i,num_data, Crhs, Wrhs, alpha, A_global, grid):
    '''Was in the file Triangle_Matrices.py. Generates the boundary h1 vector'''
    
    
    
#    A_global = A_connection(i,num_data)
#    grid = pre_processing(i, num_data)
    
    n = 2**i +1
    
    int_n = n-2
    #print int_n, 'n'
   
    
    h1 =  np.zeros([int_n**2,1])
    
    for tri in grid:
        
        
        node1 = tri[0][0]
#        ida = Coord_to_GlobalId(node1,2)

          
        node2 = tri[0][1]
#        idb = Coord_to_GlobalId(node2, 2)

              
        node3 = tri[0][2]
#        idc = Coord_to_GlobalId(node3, 2)
#        print ida, idb, idc
        
        for node in (node1,node2,node3):
            
            global_idi = Coord_to_GlobalId(node, i)
            idi = Coord_to_Id(node,i)
            #print idi, 'idi'
 
          
            
            if node[0] != 0 and node[1] != 0 and node[0] != 1 and node[1] != 1:
                #print node, [node1, node2, node3],'node'
                endpts = [node1, node2, node3]
                endpts.remove(node)
                
                for endpt in endpts:
                   
                    if endpt[0] == 0 or endpt[0] ==1 or endpt[1] == 0 or endpt[1] ==1:
                        #print endpts, 'endpts'
                        global_idj = Coord_to_GlobalId(endpt, i)
                        #print endpt, global_idj, 'endpt'
                        c = Crhs(endpt[0], endpt[1])
                        
                        w = -alpha* Wrhs(endpt[0], endpt[1])
                        #print global_idi, global_idj, 'idi, idj'
                        aentry = A_global[global_idi, global_idj]
                        
#                        print global_idi, global_idj, aentry, c
                        h1[idi] += aentry* c 
                        
                        if abs(global_idi-global_idj) == n or abs(global_idi-global_idj) ==1:
                            
                            h1[idi] += (-1)*w
                        
    return h1/2
                


def h2_bd(i,num_data, G1rhs, Wrhs, alpha, grid):
    
#    A_global = A_connection(i,num_data)
#    grid = pre_processing(i, num_data)
    
    n = 2**i +1
    
    h =1/float(n-1)
    
    int_n = n-2
    #print int_n, 'n'
   
    
    h2 =  np.zeros([int_n**2,1])
    
    for tri in grid:
        
        
        node1 = tri[0][0]
#        ida = Coord_to_GlobalId(node1,2)

          
        node2 = tri[0][1]
#        idb = Coord_to_GlobalId(node2, 2)

              
        node3 = tri[0][2]
#        idc = Coord_to_GlobalId(node3, 2)
     
        
        for node in (node1,node2,node3):
            
            global_idi = Coord_to_GlobalId(node, i)
            idi = Coord_to_Id(node,i)
            #print idi, 'idi'
 
          
            
            if node[0] != 0 and node[1] != 0 and node[0] != 1 and node[1] != 1:
                #print node, [node1, node2, node3],'node'
                endpts = [node1, node2, node3]
                endpts.remove(node)
                
                for endpt in endpts:
                   
                    if endpt[0] == 0 or endpt[0] ==1 or endpt[1] == 0 or endpt[1] ==1:
                        #print endpts, 'endpts'
                        global_idj = Coord_to_GlobalId(endpt, i)
                        #print endpt, global_idj, 'endpt'
                        g1 = G1rhs(endpt[0], endpt[1])
                        
                        w = -alpha* Wrhs(endpt[0], endpt[1])
                        #print global_idi, global_idj, 'idi, idj'
#                        aentry = A_global[global_idi, global_idj]
                        
#                        print global_idi, global_idj, aentry, w
                        if global_idi - global_idj ==n:
                            
                            h2[idi] += -2*h*w/float(6)
                            
                        if global_idi - global_idj == -n:
                            
                            h2[idi] += 2*h*w/float(6)
                            
                        if global_idi - global_idj == 1:
                            
                            h2[idi] += h*w/float(6)
                            
                        if global_idi - global_idj == -1:
                            
                            h2[idi] += -h*w/float(6)
                            
                        if global_idi - global_idj == n+1:
                            
                            h2[idi] += -h*w/float(6)
                            
                        if global_idi - global_idj ==-n-1:
                            
                            h2[idi] += h*w/float(6)
                        
                        
                        
                        if abs(global_idi-global_idj) == n or abs(global_idi-global_idj) ==1:
                            
                            h2[idi] += (-1)*alpha* g1
                        
    return h2/2
            

#h2b =  h2_bd(2,10, Xexy,XYexy, 1) 
#print h2b     
    

def h3_bd(i,num_data, G2rhs, Wrhs, alpha, grid):
    
#    A_global = A_connection(i,num_data)
#    grid = pre_processing(i, num_data)
    
    n = 2**i +1
    
    h =1/float(n-1)
    
    int_n = n-2
    #print int_n, 'n'
   
    
    h3 =  np.zeros([int_n**2,1])
    
    for tri in grid:
        
        
        node1 = tri[0][0]
        ida = Coord_to_GlobalId(node1,2)

          
        node2 = tri[0][1]
        idb = Coord_to_GlobalId(node2, 2)

              
        node3 = tri[0][2]
        idc = Coord_to_GlobalId(node3, 2)
        print ida, idb, idc
        
        for node in (node1,node2,node3):
            
            global_idi = Coord_to_GlobalId(node, i)
            idi = Coord_to_Id(node,i)
            #print idi, 'idi'
 
          
            
            if node[0] != 0 and node[1] != 0 and node[0] != 1 and node[1] != 1:
                #print node, [node1, node2, node3],'node'
                endpts = [node1, node2, node3]
                endpts.remove(node)
                
                for endpt in endpts:
                   
                    if endpt[0] == 0 or endpt[0] ==1 or endpt[1] == 0 or endpt[1] ==1:
                        #print endpts, 'endpts'
                        global_idj = Coord_to_GlobalId(endpt, i)
                        #print endpt, global_idj, 'endpt'
                        g2 = G2rhs(endpt[0], endpt[1])
                        
                        w = -alpha* Wrhs(endpt[0], endpt[1])
                        #print global_idi, global_idj, 'idi, idj'
#                        aentry = A_global[global_idi, global_idj]
                        
#                        print global_idi, global_idj, aentry, w
                        if global_idi - global_idj ==n:
                            
                            h3[idi] += 1*h*w/float(6)
                            
                        if global_idi - global_idj == -n:
                            
                            h3[idi] += -1*h*w/float(6)
                            
                        if global_idi - global_idj == 1:
                            
                            h3[idi] += -2*h*w/float(6)
                            
                        if global_idi - global_idj == -1:
                            
                            h3[idi] += 2*h*w/float(6)
                            
                        if global_idi - global_idj == n+1:
                            
                            h3[idi] += -h*w/float(6)
                            
                        if global_idi - global_idj ==-n-1:
                            
                            h3[idi] += h*w/float(6)
                        
                        
                        
                        if abs(global_idi-global_idj) == n or abs(global_idi-global_idj) ==1:
                            
                            h3[idi] += (-1)*alpha* g2
                        
    return h3/2
        
#h3b =  h3_bd(2,10, Xexy,XYexy, 1) 
#print h3b     
    

def h4_bd(i,num_data,Crhs, G1rhs, G2rhs, alpha, A_global, grid):
    
#    A_global = A_connection(i,num_data)
#    grid = pre_processing(i, num_data)
    
    n = 2**i +1
    
    h =1/float(n-1)
    
    int_n = n-2

   
    
    h4 =  np.zeros([int_n**2,1])
    
    for tri in grid:
        
        
        node1 = tri[0][0]
        ida = Coord_to_GlobalId(node1,2)

          
        node2 = tri[0][1]
        idb = Coord_to_GlobalId(node2, 2)

              
        node3 = tri[0][2]
        idc = Coord_to_GlobalId(node3, 2)
        print ida, idb, idc
        
        for node in (node1,node2,node3):
            
            global_idi = Coord_to_GlobalId(node, i)
            idi = Coord_to_Id(node,i)
            #print idi, 'idi'
 
          
            
            if node[0] != 0 and node[1] != 0 and node[0] != 1 and node[1] != 1:
                #print node, [node1, node2, node3],'node'
                endpts = [node1, node2, node3]
                endpts.remove(node)
                
                for endpt in endpts:
                   
                    if endpt[0] == 0 or endpt[0] ==1 or endpt[1] == 0 or endpt[1] ==1:
                        #print endpts, 'endpts'
                        global_idj = Coord_to_GlobalId(endpt, i)
                        #print endpt, global_idj, 'endpt'
                        c = Crhs(endpt[0], endpt[1])
                        g1 = G1rhs(endpt[0], endpt[1])
                        g2 = G2rhs(endpt[0], endpt[1])
                        
                        
                        #print global_idi, global_idj, 'idi, idj'
#                        aentry = A_global[global_idi, global_idj]
                        
#                        print global_idi, global_idj, aentry
                        ###################################################
                        # -G1
                        ####################################################
                        if global_idi - global_idj ==n:
                            
                            h4[idi] += 2*h*g1/float(6)
                            
                        if global_idi - global_idj == -n:
                            
                            h4[idi] += -2*h*g1/float(6)
                            
                        if global_idi - global_idj == 1:
                            
                            h4[idi] += -h*g1/float(6)
                            
                        if global_idi - global_idj == -1:
                            
                            h4[idi] += h*g1/float(6)
                            
                        if global_idi - global_idj == n+1:
                            
                            h4[idi] += h*g1/float(6)
                            
                        if global_idi - global_idj ==-n-1:
                            
                            h4[idi] += -h*g1/float(6)
                            
                            
                            
                        ###################################################
                        # -G2
                        ####################################################
                        if global_idi - global_idj ==n:
                            
                            h4[idi] += -1*h*g2/float(6)
                            
                        if global_idi - global_idj == -n:
                            
                            h4[idi] += 1*h*g2/float(6)
                            
                        if global_idi - global_idj == 1:
                            
                            h4[idi] += 2*h*g2/float(6)
                            
                        if global_idi - global_idj == -1:
                            
                            h4[idi] += -2*h*g2/float(6)
                            
                        if global_idi - global_idj == n+1:
                            
                            h4[idi] += h*g2/float(6)
                            
                        if global_idi - global_idj ==-n-1:
                            
                            h4[idi] += -h*g2/float(6)
                        
                        
                        
                        if abs(global_idi-global_idj) == n or abs(global_idi-global_idj) ==1:
                            
                            h4[idi] += (-1)*alpha* c
                        
    return h4/2


def Av(v):
    
    """ was in Uniform_Linearoperator.py. Not suitable for using A matrix as stencil.
    one dimensional vector -> two dimensional -> one dimensional """
    # number of interior nodes
    size = len(v)
    
    #print size
    
    # length of the interior grid
    length = int(np.sqrt(size))
    
    # Convert to two dimensional grid
    dim2_v = np.reshape(v, (length, length))
    
    # length of the entir grid
    outer_length = length+2
    
    h = 1/float(outer_length-1)
    
    # Set a dumpy whole grid
    dim2_u = np.zeros((outer_length, outer_length))
    
    dim2_u[1:-1,1:-1] = dim2_v
    
    new_v = np.zeros((outer_length, outer_length))
    
    for i in range(1, new_v.shape[0]-1):
        
        for j in range(1, new_v.shape[0]-1):
            
            new_v[i,j] = (6*dim2_u[i,j] + dim2_u[i,j+1]+ dim2_u[i,j-1] + dim2_u[i-1,j] + dim2_u[i+1, j] + dim2_u[i+1,j+1] + dim2_u[i-1,j-1])*(h**2)/float(12)
            
    assert new_v[0,0]==0, \
            "A stencil did not work on node table properly"
            
    new_v = np.reshape(new_v[1:-1,1:-1], (length**2,1))
    
    
    
    return new_v
