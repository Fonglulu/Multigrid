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
from scipy.sparse import  lil_matrix, bmat
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
    
    
def Preconditioner(grid, Ama, i, num_data,alpha):
    
#    grid = pre_processing(i, num_data)
    n = 2**i +1
    h =1/float(n-1)
    int_n = n-2
    size_L = int_n**2
#    Ama = Amatrix(i,num_data)[0]
    preG1 = lil_matrix((int_n**2, int_n**2))
    preG2 = lil_matrix((int_n**2, int_n**2))
    Lma = lil_matrix((int_n**2, int_n**2))
    
    count = 0
    for tri in grid:

               count+=1
               #print count
               # node 1 should be the right angle node
               node1 = tri[0][0]
               idi = Coord_to_Id(node1,i)
               
               node2 = tri[0][1]
               idj = Coord_to_Id(node2,i)
              
               node3 = tri[0][2]       
               idk = Coord_to_Id(node3,i)
#               print (node1,node2, node3), (idi, idj, idk), 'tri'
               
               if 0<=idi<= size_L and 0<= idj<=size_L and 0<=idk<=size_L:
#                   print (node1,node2, node3), (idi, idj, idk), 'tri'
                   
                   Lma[idi, idi] =4
                   Lma[idj, idj] =4
                   Lma[idk, idk] =4
              
   
               
                   if abs(idi-idj) ==1:
                    
                       
                       preG2[idi, idj] =1*h
                       preG2[idj, idi] = -1*h
                       
                       Lma[idi, idj] =-1
                       Lma[idj, idi] = -1 
                   

                       
                   if abs(idj-idk) ==1 :
                       
                       preG2[idj,idk] = 1*h
                       preG2[idk,idj] =  -1*h
                       
                       Lma[idj,idk] =-1
                       Lma[idk,idj] =-1
#                       
                   if abs(idi-idj) == int_n:
               
                       preG1[idi, idj] = 1*h
                       preG1[idj, idi] = -1*h
                       
                       Lma[idi, idj] =-1
                       Lma[idj,idi] =-1
                       
                   if abs(idj-idk) == int_n:
                       
                        preG1[idj,idk] =  1*h
                        preG1[idk, idj] = -1*h
                        
                        Lma[idj,idk] = -1
                        Lma[idk, idj] = -1
                        
    S_sym = bmat([[Ama, None, None,  np.sqrt(alpha)*Lma],\
                        [None, Lma, None, -preG1.T],\
                        [None, None, Lma,  -preG2.T],\
                        [ np.sqrt(alpha)*Lma, -preG1, -preG2, None]])


               
    return S_sym




    
    
    

def Amatrix(grid, i, num_data):
 
#    grid = pre_processing(i, num_data)

    
    n = 2**i +1
    
    int_n = n-2
    #print int_n, 'n'
    size_A = int_n**2

  
    Amatrix = lil_matrix((int_n**2, int_n**2))
    
    Dvector =  np.zeros((int_n**2))


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
                       
                        Dvector[idi] += i_data
                        
                    if (0<=idj<=size_A): 
                       
                        Amatrix[idj, idj] +=jj
                        
                        Dvector[idj] += j_data
          
                    if (0<=idk<=size_A):
                        
                        Amatrix[idk, idk] += kk
                    
                        Dvector[idk] += k_data
        
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


def A_connection(grid, i, num_data):
    
#    start = time.time()
#    grid = pre_processing(i, num_data)
#    done = time.time()
#    elapsed = done - start
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
                   
                  
              
                    jj = Polynomial_eval(node2, node1, node3, data) * Polynomial_eval(node2, node1, node3, data)

            
                    kk = Polynomial_eval(node3, node1, node2, data) * Polynomial_eval(node3, node1, node2, data)

                
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







    
    
def h_boundary(grid, A_global, d, i,num_data, Crhs,G1rhs, G2rhs, Wrhs, alpha):
    
    
#    grid = pre_processing(i, num_data)
#    
#    A_global = A_connection(grid, i,num_data)
#    
#    d = Amatrix(grid, i, num_data)[-1]
    
    n = 2**i +1
    
    int_n = n-2
    
    rhs =  np.zeros((4*int_n**2,))

    h =1/float(n-1)
    
    h1 =  np.zeros((int_n**2,))
    h2 =  np.zeros((int_n**2,))
    h3 =  np.zeros((int_n**2,))
    h4 =  np.zeros((int_n**2,))
    
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
                            
                            h4[idi] += (-1)* c
                            
    rhs[0:int_n**2] = d - h1/2
    rhs[int_n**2:2*int_n**2] = -h2/(2*float(np.sqrt(alpha)))
    rhs[2*int_n**2:3*int_n**2]=-h3/(2*float(np.sqrt(alpha)))
    rhs[3*int_n**2:] =-h4*float(np.sqrt(alpha))/2
    
    
    return rhs





    