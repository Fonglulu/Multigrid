#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:04:57 2019

@author: shilu
"""

from grid.NodeTable import node_iterator, not_slave_node
from scipy.sparse import csc_matrix, lil_matrix, bmat, coo_matrix, csr_matrix
from grid.ConnectTable import connect_iterator 
import numpy as np
import math
from BuildEquation import build_equation_linear_2D, set_polynomial_linear_2D,\
Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY, NIn_triangle, build_matrix_fem_2D
from Triangle import triangle_iterator, Interior_triangle_iterator

def Amatrix(grid, Coord, Nl, h):
    
            
    Amatrix = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
        
        if not node.get_slave():
            
            idi = int(node.get_value())
            
          
            
            #print idi, node.get_coord()

            
            #print node.get_value(), node.get_node_id()._id_no
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
                idj = int(grid.get_value(endpt)) 
               
                
                #print idj, grid.get_coord(endpt)
                
                #print endpt._id_no, grid.get_value(endpt)
                
                aentry = grid.get_matrix_value(node.get_node_id(), endpt)[1]

                #print aentry
                
                if not grid.get_slave(endpt):
                    
                    Amatrix[idi, idj] =aentry


    return Amatrix/float(len(Coord))

def Stencil_A(grid, Nl, h):
       
    Amatrix = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
        
        if not node.get_slave():
            
            idi = int(node.get_value())
            
            
            
            # Loop over the matrix entries in the current row
            for endpt in connect_iterator(grid, node.get_node_id()):
                
                # Which column corresponds to the current node?
                idj = int(grid.get_value(endpt))
                
                
                # We must not include slave nodes in the matrix columns
                if not grid.get_slave(endpt):
                    
                    if idi - idj == 0:
                        
                        Amatrix[idi, idj]= h**2/2
                        
                    elif abs(idi - idj) == 1:
                        
                        Amatrix[idi, idj] = h**2/12
                    
                    elif abs(idi - idj) == int(np.sqrt(len(Nl))):
                        
                        Amatrix[idi, idj] = h**2/12
                        
                    elif abs(idi - idj) == int(np.sqrt(len(Nl)))+1:
                        
                        Amatrix[idi, idj] = h**2/12
                        
    return Amatrix
                        
                    
                    

def Lmatrix(grid, Nl):
    
    Lmatrix = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt1)[0] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        Lmatrix[i, j] = lentry
    return Lmatrix



def Stencil_L(grid, Nl):
    
    
    Lmatrix = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            idi = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    idj = int(grid.get_value(endpt1))

                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        
                        if idi - idj == 0 :
                            
                            Lmatrix[idi, idj] = 4
                                
                        elif abs(idi-idj) == 1:
                            
                            Lmatrix[idi, idj] = -1
                            
                        elif abs(idi - idj) == int(np.sqrt(len(Nl))):
                            
                            Lmatrix[idi, idj] = -1
                            
                            
                            
                            
                            
                            
    return Lmatrix
    


#Generate G1 matrix
def G1(grid, Nl):
    


    
    G1 = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
            #print node.get_coord(), i
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    #print  j
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    g1 = grid.get_matrix_value(node.get_node_id(), endpt1)[2] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        G1[i, j] = g1
    return G1



def FD_G1(grid, Nl, h):
    
    
    G1 = csr_matrix((len(Nl),len(Nl)))
    
#    for node in Nl:
#        
#        node.set_value(Nl.index(node))
        
    for node in node_iterator(grid):
        
        if not node.get_slave():
            
            i = int(node.get_value())
            
            for endpt1 in connect_iterator(grid, node.get_node_id()):
                
                j = int(grid.get_value(endpt1))
                
                if not grid.get_slave(endpt1):
                    
                    if i-j == 0:
                        
                        G1[i,j] = 2*h
                        
                    elif i -j == np.sqrt(len(Nl)):
                        
                        G1[i,j] = -h
                        
                    elif i -j == -np.sqrt(len(Nl)):
                        
                        G1[i,j] = -h
                        
#                    elif i-j == 1:
#                        
#                        
#                        G1[i,j] = -h
#                        
#                    elif i-j == -1:
#                        
#                        
#                        G1[i,j] = -h
#                        
#                    elif i-j == np.sqrt(len(Nl))+1:
#                        
#                        G1[i,j] =-h
#                        
#                    elif i-j == -np.sqrt(len(Nl))+1:
#                        
#                        G1[i,j] = h
        
        
    return G1



def G2(grid, Nl):
    
    

    G2 = csr_matrix((len(Nl), len(Nl)))
    
    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    g2 = grid.get_matrix_value(node.get_node_id(), endpt1)[3] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        G2[i, j] = g2
                    
                    
    return G2


def FD_G2(grid, Nl, h):
    

    
    G2 = csr_matrix((len(Nl), len(Nl)))
    
    
    
#    for node in Nl:
#        
#        node.set_value(Nl.index(node))
        
    for node in node_iterator(grid):
        
        if not node.get_slave():
            
            i = int(node.get_value())
            
            for endpt1 in connect_iterator(grid, node.get_node_id()):
                
                j = int(grid.get_value(endpt1))
                
                if not grid.get_slave(endpt1):
                    
                    if i-j == 0:
                        
                        G2[i,j] = 2*h
                        
                    elif abs(i -j) == 1:
                        
                        G2[i,j] = -h
                        
                        
#                    elif abs(i -j) == np.sqrt(len(Nl)):
#                        
#                        G2[i,j] = -h
#                    
#                    elif abs(i-j) == np.sqrt(len(Nl))+1:
#                        
#                        G2[i,j] = -h

        
        
    return G2

def dvector(grid, data, Nl, Coord):
    


               
    dvector = np.zeros([len(Nl),1])
    
    for tri in Interior_triangle_iterator(grid):
        
        if tri[1].get_node_id() < tri[2].get_node_id():
            
            basis1 = set_polynomial_linear_2D(tri[0], tri[2], tri[1])
            
            Idi = int(tri[0].get_value())
            
            
            for i in range(len(Coord)):
                
                
                if NIn_triangle(tri[0] , tri[1], tri[2], Coord[i]):
                    
                        
                        
                        dvector[Idi,0] += basis1.eval(Coord[i][0], Coord[i][1]) * data[i]
                        
                        
                    
    return dvector/float(len(Coord))

def h1_bd(grid, Crhs, wrhs, Nl, Coord): 
    
    
    
    h1 = np.zeros((len(Nl), 1))
        
    
    for node in not_slave_node(grid):
        
        # Ignore slave (or boundary) nodes
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            #print node.get_value(), node.get_node_id()._id_no
            
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
                #j = int(grid.get_value(endpt))
                
                #print j
                
            
                if grid.get_slave(endpt):
                    
                    
                    coord = grid.get_coord(endpt)
            
                    c = Crhs(coord[0], coord[1])
                    
                    w = wrhs(coord[0], coord[1])

                                        
                    #print aentry
                    
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
            
                    aentry = grid.get_matrix_value(node.get_node_id(), endpt)[1]/float(len(Coord))
            
                    h1[i] += c* aentry + w * lentry
    return h1


def h2_bd(grid, g1rhs, wrhs, alpha, Nl):
    

    
    
    h2= np.zeros([len(Nl), 1])
    
    for node in not_slave_node(grid):
        
        # Ignore slave (or boundary) nodes
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
            
                if grid.get_slave(endpt):
                    
                    coord = grid.get_coord(endpt)
                    
                    #print coord
            
                    g1 = g1rhs(coord[0], coord[1])
                    #print g1
                    
                    w = wrhs(coord[0], coord[1])

                    
                    #print endpt._id_no
            
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
                    
                    g1entry = grid.get_matrix_value(node.get_node_id(), endpt)[2]
                    
                    #print aentry
            
                    #h2[i] += alpha * g1 * lentry - w * g1entry #G1, G2
                   
                    h2[i] += alpha * g1 * lentry +w * g1entry  #-G1, -G2
                    
    return h2


def h3_bd(grid, g2rhs, wrhs, alpha, Nl):
    

    
    
    h3= np.zeros([len(Nl), 1])
    
    for node in not_slave_node(grid):
        
        # Ignore slave (or boundary) nodes
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
            
                if grid.get_slave(endpt):
                    
                    coord = grid.get_coord(endpt)
            
                    g2 = g2rhs(coord[0], coord[1])
                    
                    w = wrhs(coord[0], coord[1])
                    
                    #print c, w
                    
                    #print endpt._id_no
            
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
                    
                    # G2 entries on boundary
                    g2entry = grid.get_matrix_value(node.get_node_id(), endpt)[3]
                    
                    #print aentry
                    
                    #h3[i] += alpha * g2* lentry - w* g2entry # G1, G2
                    h3[i] += alpha * g2* lentry + w* g2entry  # -G1, -G2
                    
    return h3


def h4_bd(grid, Crhs, g1rhs, g2rhs, Nl):
    

    
    
    h4 = np.zeros([len(Nl),1])
    
    for node in not_slave_node(grid):
    
    
            # Which row corresponds to the current node?
            i = int(node.get_value())
            
            for endpt in connect_iterator(grid, node.get_node_id()):
                
            
                if grid.get_slave(endpt):
                    
                    coord = grid.get_coord(endpt)
                    
            
                    c = Crhs(coord[0], coord[1])
                    
                    g1 = g1rhs(coord[0], coord[1])
                    
                    g2 = g2rhs(coord[0], coord[1])
                    
                    #print c
                    
                    #print endpt._id_no
            
                    lentry = grid.get_matrix_value(node.get_node_id(), endpt)[0]
                    
                    g1entry = grid.get_matrix_value(node.get_node_id(), endpt)[2]
                    
                    g2entry = grid.get_matrix_value(node.get_node_id(), endpt)[3]
                    
                    #print aentry
                    
                    
            
                    #h4[i] += c* lentry + g1entry * g1  + g2entry * g2 # G1, G2
                    h4[i] += c* lentry - g1entry * g1  - g2entry * g2 # -G1, -G2
                    
                    #print h4
                    
    return h4