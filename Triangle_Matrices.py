#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:05:45 2019

@author: shilu
"""

#@profile   

import time
import numpy as np 
from PreProcess import pre_processing, Polynomial_eval
from scipy.sparse import  lil_matrix
from functions import Exy, XYexy, Xexy, Yexy


def Coord_to_Id(coord,i):
    
    n = 2**i+1
    
    int_n = n-2
    
    h =1/float(n-1)
    
    
    
    num_x = int(coord[0]/float(h))

    assert abs(num_x - coord[0]/float(h)) <= 1e-8, 'coordinate is not in the grid'
    
    num_y = int(coord[1]/float(h))
    


    if (num_x != 0) and (num_y != n-1) and (num_x != n-1) and (num_y != 0):
                
                int_id  = (num_x-1) *int_n + (num_y-1)
                
    else:
                
                int_id =-1
    
    return int_id



def Coord_to_GlobalId(coord,i):
    
    n = 2**i+1

    
    h =1/float(n-1)
    
    num_x = int(coord[0]/float(h))
    
    num_y = int(coord[1]/float(h))
    
    global_id = num_x*n+num_y
    
    return global_id
    
    
    
    

def Amatrix(i, num_data):
    start = time.time()
    grid = pre_processing(i, num_data)
    done = time.time()
    elapsed = done - start
    print(elapsed)
    
    n = 2**i +1
    
    int_n = n-2
    #print int_n, 'n'
    size_A = int_n**2

  
    Amatrix = lil_matrix((int_n**2, int_n**2))
    
    Dvector =  np.zeros([int_n**2,1])


    for tri in grid:



               # node 1 should be the right angle node
               node1 = tri[0][0]
               idi = Coord_to_Id(node1,i)
          
               node2 = tri[0][1]
               idj = Coord_to_Id(node2,i)
              
               node3 = tri[0][2]       
               idk = Coord_to_Id(node3,i)
               
              

               for data in tri[1]._data:
                   
 
#                    print 'evaluation'
#                    start = time.time()
                    ii = Polynomial_eval(node1, node2, node3, data) * Polynomial_eval(node1, node2, node3, data)
                   
                    i_data = Polynomial_eval(node1, node2, node3, data) * data[2]
              
                    jj = Polynomial_eval(node2, node1, node3, data) * Polynomial_eval(node2, node1, node3, data)
                    
                    j_data = Polynomial_eval(node2, node1, node3, data) * data[2]
            
                    kk = Polynomial_eval(node3, node1, node2, data) * Polynomial_eval(node3, node1, node2, data)
                    
                    k_data = Polynomial_eval(node3, node1, node2, data) * data[2]
                
                    ij = Polynomial_eval(node1, node2, node3, data) * Polynomial_eval(node2, node1, node3, data)
            
                    ik = Polynomial_eval(node1, node2, node3, data) * Polynomial_eval(node3, node1, node2, data)
                
                    jk =  Polynomial_eval(node2, node1, node3, data) * Polynomial_eval(node3, node1, node2, data)
#                    done = time.time()
#                    elapsed = done - start
#                    print(elapsed)
#                    print 'evaluation'
#                    print 'allocation'
#                    start = time.time()
                   
                    if (0<=idi <=size_A) :
                     
                        Amatrix[idi,idi] += ii
                       
                        Dvector[idi,0] += i_data
                        
                    if (0<=idj<=size_A): 
                       
                        Amatrix[idj, idj] +=jj
                        
                        Dvector[idj,0] += j_data
          
                    if (0<=idk<=size_A):
                        
                        Amatrix[idk, idk] += kk
                    
                        Dvector[idk,0] += k_data
        
                    if (0<=idi <=size_A) and (0<=idj<= size_A):
                        
                        
                        Amatrix[idi, idj] += ij
                
                        Amatrix[idj, idi]+= ij
                        
                    if (0<=idi<=size_A) and (0<=idk<=size_A):
                        
                        
               
                        Amatrix[idi,idk] += ik
                    
                        Amatrix[idk,idi] += ik
                    if (0<=idj <=size_A) and (0<=idk <=size_A):
                        
                        
                        
                        Amatrix[idj, idk] += jk
                
                        Amatrix[idk, idj] += jk
                    
#                    done = time.time()
    
#                    elapsed = done - start
#                    print(elapsed)
                    
                  
    return Amatrix/(num_data**2), Dvector/(num_data**2)


def A_connection(i, num_data):
    
    start = time.time()
    grid = pre_processing(i, num_data)
    done = time.time()
    elapsed = done - start
    #print(elapsed)
    
    n = 2**i +1

    #print int_n, 'n'
    size_A = n**2

  
    Amatrix = lil_matrix((n**2, n**2))
    


    for tri in grid:



               # node 1 should be the right angle node
               node1 = tri[0][0]
               idi = Coord_to_GlobalId(node1,i)
          
               node2 = tri[0][1]
               idj = Coord_to_GlobalId(node2,i)
              
               node3 = tri[0][2]       
               idk = Coord_to_GlobalId(node3,i)
               
              

               for data in tri[1]._data:
                   
 
#                    print 'evaluation'
#                    start = time.time()
                    ii = Polynomial_eval(node1, node2, node3, data) * Polynomial_eval(node1, node2, node3, data)
                   
                    i_data = Polynomial_eval(node1, node2, node3, data) * data[2]
              
                    jj = Polynomial_eval(node2, node1, node3, data) * Polynomial_eval(node2, node1, node3, data)
                    
                    j_data = Polynomial_eval(node2, node1, node3, data) * data[2]
            
                    kk = Polynomial_eval(node3, node1, node2, data) * Polynomial_eval(node3, node1, node2, data)
                    
                    k_data = Polynomial_eval(node3, node1, node2, data) * data[2]
                
                    ij = Polynomial_eval(node1, node2, node3, data) * Polynomial_eval(node2, node1, node3, data)
            
                    ik = Polynomial_eval(node1, node2, node3, data) * Polynomial_eval(node3, node1, node2, data)
                
                    jk =  Polynomial_eval(node2, node1, node3, data) * Polynomial_eval(node3, node1, node2, data)
#                    done = time.time()
#                    elapsed = done - start
#                    print(elapsed)
#                    print 'evaluation'
#                    print 'allocation'
#                    start = time.time()
                   
                    if (0<=idi <=size_A) :
                
                        Amatrix[idi,idi] += ii
                       
                        #Dvector[idi,0] += i_data
                        
                    if (0<=idj<=size_A): 
                       
                        Amatrix[idj, idj] +=jj
                        
                        #Dvector[idj,0] += j_data
          
                    if (0<=idk<=size_A):
                        
                        Amatrix[idk, idk] += kk
                    
                        #Dvector[idk,0] += k_data
        
                    if (0<=idi <=size_A) and (0<=idj<= size_A):
                        
                        
                        Amatrix[idi, idj] += ij
                
                        Amatrix[idj, idi]+= ij
                        
                    if (0<=idi<=size_A) and (0<=idk<=size_A):
                        
                        
               
                        Amatrix[idi,idk] += ik
                    
                        Amatrix[idk,idi] += ik
                    if (0<=idj <=size_A) and (0<=idk <=size_A):
                        
                
                        Amatrix[idj, idk] += jk
                
                        Amatrix[idk, idj] += jk
                    
#                    done = time.time()
    
#                    elapsed = done - start
#                    print(elapsed)
                    
                  
    return Amatrix.todense()/(num_data**2)


def h1_bd(i,num_data, Crhs, Wrhs, alpha, A_global, grid):
    
    
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

#h4b = h4_bd(2, 10, Exy, Xexy, Yexy,1)
#print h4b



def rhs(i, num_data, alpha, Crhs, G1rhs, G2rhs, Wrhs):
    
    A_global = A_connection(i,num_data)
    
    grid =  pre_processing(i, num_data)
    
    n = 2**i +1
    
#    h =1/float(n-1)
    
    int_n = n-2
    
    rhs =  np.zeros([4*int_n**2,1])
    
    h1 = h1_bd(i,num_data, Crhs, Wrhs, alpha,A_global, grid)
    
    h2 = h2_bd(i,num_data,G1rhs,Wrhs,alpha,grid)
    
    h3 = h3_bd(i, num_data, G2rhs,Wrhs, alpha, grid)
    
    h4 = h4_bd(i,num_data, Crhs, G1rhs, G2rhs, alpha, A_global, grid)
    
    rhs[0:int_n**2] = h1
    rhs[int_n**2:2*int_n**2] = h2
    rhs[2*int_n**2:3*int_n**2]=h3
    rhs[3*int_n**2:] =h4
    
    return rhs

rhsb = rhs(2, 10, 1, Exy, Xexy, Yexy, XYexy)

print rhsb    
    
    
def h_boundary(i,num_data, Crhs,G1rhs, G2rhs, Wrhs, alpha):
    
    
    A_global = A_connection(i,num_data)
    grid = pre_processing(i, num_data)
    
    n = 2**i +1
    
    int_n = n-2

    h =1/float(n-1)
    
    h1 =  np.zeros([int_n**2,1])
    h2 =  np.zeros([int_n**2,1])
    h3 =  np.zeros([int_n**2,1])
    h4 =  np.zeros([int_n**2,1])
    
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

                        global_idj = Coord_to_GlobalId(endpt, i)

                        c = Crhs(endpt[0], endpt[1])
                        
                        g1 = G1rhs(endpt[0], endpt[1]) 
                
                        g2 =  G2rhs(endpt[0], endpt[1])
                        
                        w = -alpha* Wrhs(endpt[0], endpt[1])
               
                        aentry = A_global[global_idi, global_idj]
                        
                        ###################################################
                        # A  
                        ####################################################  

                        h1[idi] += aentry* c 
                                                    
                            
                        ###################################################
                        # -G1^T = G1
                        ####################################################  
                            
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
                        ###################################################
                        # -G2^T = G2
                        ####################################################                            
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
                            
                            
                            
                        ##################################################
                        # -G2
                        ###################################################
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
                            

                        
                        ##################################################
                        # L
                        ###################################################
                        
                        
                        
                        
                        if abs(global_idi-global_idj) == n or abs(global_idi-global_idj) ==1:
                            
                            
                            h1[idi] += (-1)*w 
                            
                            h2[idi] += (-1)*alpha* g1
                            
                            h3[idi] += (-1)*alpha* g2
                            
                            h4[idi] += (-1)*alpha* c
    return h1/2, h2/2, h3/2, h4/2


    