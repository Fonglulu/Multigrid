#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:50:19 2019

@author: shilu
"""
import time 
import numpy as np
from functions import Linear, Linear2x
np.set_printoptions(precision=4)


i = 2

n= 2**i+1

# Find the spacing
h=1/float(n-1)

# Set the mesh grid, that is interior
#x1, y1 = np.meshgrid(np.arange(h, 1, h), np.arange(h, 1, h))
x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))

nodes = np.vstack([x1.ravel(), y1.ravel()]).T



# Set up data
x = np.linspace(0, 1.0,100)
y = np.linspace(0, 1.0,100)
X, Y = np.meshgrid(x,y)

data = Linear2x(X,Y)
data = data.flatten()

coordx = X.flatten()
coordy = Y.flatten()
Coord = zip(coordx, coordy)


def Polynomial_eval(node1, node2, node3, data_coord):
    
    from math import fabs
    #print data_coord
    x1 = node1[0]
    y1 = node1[1]
    x2 = node2[0]
    y2 = node2[1]
    x3 = node3[0]
    y3 = node3[1]
    
    division = (y1-y2)*(x2-x3)-(y2-y3)*(x1-x2);
    #print division
    
    
    assert fabs(division) > 1.0E-12, "divide by zero"
    
    const = (x3*y2 - y3*x2)/division
    xcoe = (y3-y2)/division
    ycoe = (x2-x3)/division
    
    return data_coord[0]*xcoe+data_coord[1]*ycoe+const


def In_triangle(node1, node2, node3, data_coord):
    
    value1 = Polynomial_eval(node1, node2, node3, data_coord)
    value2 = Polynomial_eval(node2, node1, node3, data_coord)
    value3 = Polynomial_eval(node3, node2, node1, data_coord)
    
    return  (value1 >=0.0 and value2>=0.0 and value3 >=0.0)


def  dvector(Coord, data, nodes,n):
    
    np.set_printoptions(precision=16)
    from scipy import zeros
    from copy import copy
    # initilise the dvector.
    qdvector = zeros([(n)**2,1])
    
    for i in range(len(Coord)):
        
    
        
        data_i = data[i]
        data_coord = Coord[i]
        #print Coord[i]
        
        # difference consists the distances between give data and all grid nodes
        
        difference = Coord[i] - nodes
        
        
        
        distance = difference[:,0]**2 +difference[:,1]**2
        
        difference = difference.tolist()
        
        distance = distance.tolist()
        
        #print distance, 'distance'
        
        
        node_list = copy(nodes.tolist())
        
        #print node_list, 'node_list'
        
        #print node_list
        
        
        # Find the closest node
        #  the index of  smallest distance 
        node1 = node_list[distance.index(sorted(distance)[0])]
       
        #print distance.index(sorted(distance)[0]), node1, 'node1'
        
        
        
        # The index of cloest node
        IDi = node_list.index(node1)
        #print IDi, 'IDi'
        
        # Find the second & third closest node
        
        # Find the index of second smallest distance
        node2 = node_list[distance.index(sorted(distance)[1])]
       
        
        index2 = distance.index(sorted(distance)[1])
        #print distance.index(sorted(distance)[1]), node2, 'node2'
        
        IDj = node_list.index(node2)
        #print IDj, 'IDj'
        
        # In case the seond and the third node have the same distance to the coordinate
        # And the second one had been picked up again to be the third.
        distance[index2] = distance[index2]+3
        
        
        
        #distance.pop(distance.index(sorted(distance)[1]))
        #print distance
      
        
        # Find the index of third smallest distance
        #node_list_2 = copy(node_list)
        #node_list_2.remove(node2)
        node3 = node_list[distance.index(sorted(distance)[1])]
        #print distance.index(sorted(distance)[1]), node3, 'node3'
        
        
        IDk = node_list.index(node3)
        #print IDk, 'IDk'
        
        #print node1, node2, node3, Coord[i]
        
        if In_triangle(node1, node2, node3, Coord[i]):
          
        
            qdvector[IDi,0] += Polynomial_eval(node1, node2, node3, data_coord)* data_i
            
            qdvector[IDj,0] += Polynomial_eval(node2, node1, node3, data_coord)* data_i
            
            qdvector[IDk,0] += Polynomial_eval(node3, node1, node2, data_coord)* data_i
        
        
    return qdvector

#start = time.time()
fast_rhs =dvector(Coord, data, nodes, n)/float(len(Coord))
#
fast_rhs = np.reshape(fast_rhs, (n,n))[1:-1,1:-1]
fast_rhs= np.reshape(fast_rhs, ((n-2)**2,1))
#done = time.time()
#elapsed = done - start
#print elapsed
#        
        
       
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    