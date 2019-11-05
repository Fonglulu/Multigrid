#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:42:18 2019

@author: shilu
"""

import numpy as np
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import csr_matrix, bmat
from scipy import eye






def G1v(v):
    
    size = len(v)
    
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
    
    
    for x in range(1, new_v.shape[0]-1):
        
        for y in range(1, new_v.shape[0]-1):
            
            # My version
            new_v[x,y] = (-2*dim2_u[x-1, y]+2*dim2_u[x+1,y] - dim2_u[x, y+1] + dim2_u[x, y-1] +dim2_u[x+1, y+1] - dim2_u[x-1, y-1])*h/float(6)
            
            
    assert new_v[0,0]==0, \
            "G1 stencil did not work on node table properly"
            
    new_v = np.reshape(new_v[1:-1,1:-1], (length**2,1))
    
    
    
    return new_v


def G1Tv(v):
    
    size = len(v)
    
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
    
    
    for x in range(1, new_v.shape[0]-1):
        
        for y in range(1, new_v.shape[0]-1):
            
            # My version
            new_v[x,y] = (2*dim2_u[x-1,y]-2*dim2_u[x+1,y] + dim2_u[x, y+1] - dim2_u[x, y-1] -dim2_u[x+1, y+1] + dim2_u[x-1, y-1])*h/float(6)
            
            
    assert new_v[0,0]==0, \
            "G1 stencil did not work on node table properly"
            
    new_v = np.reshape(new_v[1:-1,1:-1], (length**2,1))
    
    
    
    return new_v
    


def G2v(v):
    
    size = len(v)
    
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
    
    
    for x in range(1, new_v.shape[0]-1):
        
        for y in range(1, new_v.shape[0]-1):
            
            # My version
           new_v[x,y] = (dim2_u[x-1,y] - dim2_u[x+1,y] +2 * dim2_u[x,y+1] -2*dim2_u[x,y-1] - dim2_u[x-1, y-1] +dim2_u[x+1,y+1])*h/float(6)
            
            
    assert new_v[0,0]==0, \
            "G2 stencil did not work on node table properly"
            
    new_v = np.reshape(new_v[1:-1,1:-1], (length**2,1))
    
    
    
    return new_v



def G2Tv(v):
    
    size = len(v)
    
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
    
    
    for x in range(1, new_v.shape[0]-1):
        
        for y in range(1, new_v.shape[0]-1):
            
            # My version
           new_v[x,y] = (-dim2_u[x-1,y] + dim2_u[x+1,y] - 2 * dim2_u[x,y+1] + 2*dim2_u[x,y-1] + dim2_u[x-1, y-1] - dim2_u[x+1,y+1])*h/float(6)
            
            
    assert new_v[0,0]==0, \
            "G2 stencil did not work on node table properly"
            
    new_v = np.reshape(new_v[1:-1,1:-1], (length**2,1))
    
    
    
    return new_v
    


def Lv(v):
    
    size = len(v)
    
    # length of the interior grid
    length = int(np.sqrt(size))
    
    # Convert to two dimensional grid
    dim2_v = np.reshape(v, (length, length))
    
    # length of the entir grid
    outer_length = length+2

    
    # Set a dumpy whole grid
    dim2_u = np.zeros((outer_length, outer_length))
    
    dim2_u[1:-1,1:-1] = dim2_v
    
    new_v = np.zeros((outer_length, outer_length))
    
    
    for i in range(1, new_v.shape[0]-1):
        
        for j in range(1, new_v.shape[0]-1):
            
            # My version
           new_v[i,j] = (4* dim2_u[i,j] -  dim2_u[i-1,j] - dim2_u[i+1,j] - dim2_u[i,j-1] - dim2_u[i, j+1])
            
            
    assert new_v[0,0]==0, \
            "L stencil did not work on node table properly"
            
    new_v = np.reshape(new_v[1:-1,1:-1], (length**2,1))
    
    
    
    return new_v
    
    
    
    
    
    
    
    

def Sv(v):
    """ Take v as one dimensional vector """
    
    
    size = int(len(v)/float(4))
    
    #print Av(v[0:size])+ Lv(v[3*size:])
    #print v.shape, 'v'
    new_v = np.zeros(v.shape[0])
   
    
    #length = int(np.sqrt(size))
    
    new_v[0:size] = A_LinearOperator.matvec(v[0:size]) + L_LinearOperator.matvec(v[3*size:])
    
    new_v[size: 2*size] = np.sqrt(beta)* L_LinearOperator.matvec(v[size:2*size]) + G1_LinearOperator.rmatvec(v[3*size:])
    
    new_v[2*size: 3*size] = np.sqrt(beta) * L_LinearOperator.matvec(v[2*size : 3*size]) + G2_LinearOperator.rmatvec(v[3*size:])
    
    new_v[3*size:] = L_LinearOperator.matvec(v[0:size]) +G1_LinearOperator.matvec(v[size:2*size]) + G2_LinearOperator.matvec(v[2*size:3*size])
    
    
    return new_v
##    
#  
#    
    
    
    
    


#def Setup_LinearOperators(i, alpha):
#    
#  
#    global A_LinearOperator, L_LinearOperator, G1_LinearOperator, G2_LinearOperator, beta
#    
#    beta =alpha 
#    
#    n = 2**i+1
#    
#    
#    h =1/float(n-1)
#
#    _value = np.zeros((4,n-2,n-2))
#    
#    int_c = np.zeros(((n-2)**2, (n-2)**2))
#    
#    print int_c.shape
#    
#    int_g1 = np.zeros(((n-2)**2, (n-2)**2))
#    
#    int_g2 = np.zeros(((n-2)**2, (n-2)**2))
#    
#    int_w = np.zeros(((n-2)**2, (n-2)**2))
#    
#    A_LinearOperator = LinearOperator(int_c.shape, matvec = Av, rmatvec = Av)
#    
#    L_LinearOperator = LinearOperator(int_c.shape, matvec = Lv, rmatvec = Lv)
#    
#    G1_LinearOperator = LinearOperator(int_c.shape, matvec = G1v, rmatvec = G1Tv)
#    
#    G2_LinearOperator = LinearOperator(int_c.shape, matvec = G2v, rmatvec = G2Tv)
#    
#    S_LinearOperator = LinearOperator(tuple(4*x for  x in int_c.shape), matvec = Sv, rmatvec = Sv)
#    
#    
# 
#    return S_LinearOperator



def Setup_LinearOperator(Amatrix, i,num_data, alpha):
    
#    from Triangle_Matrices import Amatrix
#    
#    A = Amatrix(i, num_data)[0]
   
    
    n = 2**i+1
    
    int_c = np.zeros(((n-2)**2, (n-2)**2))
    
    L_LinearOperator = LinearOperator(int_c.shape, matvec = Lv, rmatvec = Lv)
    
    G1_LinearOperator = LinearOperator(int_c.shape, matvec = G1v, rmatvec = G1Tv)
    
    G2_LinearOperator = LinearOperator(int_c.shape, matvec = G2v, rmatvec = G2Tv)
    
    def Sv(v):
        
  
        size = int(len(v)/float(4))
        
        new_v = np.zeros(v.shape[0])
       
        new_v[0:size] = Amatrix.dot(v[0:size]) + np.sqrt(alpha)*L_LinearOperator.matvec(v[3*size:])

        new_v[size: 2*size] =  L_LinearOperator.matvec(v[size:2*size]) + G1_LinearOperator.matvec(v[3*size:])
    
        new_v[2*size: 3*size] =  L_LinearOperator.matvec(v[2*size : 3*size]) + G2_LinearOperator.matvec(v[3*size:])
    
        new_v[3*size:] = np.sqrt(alpha)*L_LinearOperator.matvec(v[0:size]) -G1_LinearOperator.matvec(v[size:2*size]) - G2_LinearOperator.matvec(v[2*size:3*size])
        
        
        return new_v
    S_LinearOperator = LinearOperator(tuple(4*x for  x in int_c.shape), matvec =Sv, rmatvec =Sv)
    
    return S_LinearOperator

def negh(n,  idi):
    
    
    # number of nodes on edge
    int_n = n-2
    
    # Find the divisor, which tells the rows
    x_ind = idi % int_n 

    
    # Find the remainder, which tells the cols
    y_ind = idi // int_n 
    
    negh_set = []
    
    for k in [-1,0,1]:
        
        if x_ind + k in [x for x in range(0, int_n)]:
            

            for j in [-1,0,1]:
                
                if y_ind + j in [ x for x in range(0,int_n)]:
                    

                    
                    negh_id = (y_ind+j)* int_n+ x_ind+k
                    

                    
                    negh_set.append(negh_id)
                    
    return negh_set
                    
                    
            
            

def Stencil_A(i):
    
    # Total number of nodes
    
    n = 2**i+1
    
    h =1/float(n-1)
       
    #Amatrix = csr_matrix(((n-2)**2, (n-2)**2))
    Amatrix = eye((n-2)**2)
    
    for idi in range(0, (Amatrix.shape[0])):
        
        
        negh_set = negh(n,  idi)
        

                
        for idj in negh_set:
                    
                    
                    if idi - idj == 0:
                        
                        Amatrix[idi, idj]= h**2/2
                        
                    elif abs(idi - idj) == 1:
                        
                        Amatrix[idi, idj] = h**2/12
                    
                    elif abs(idi - idj) == int(n-2):
                        
                        Amatrix[idi, idj] = h**2/12
                        
                    elif abs(idi - idj) == int(n-2)+1:
                        
                        Amatrix[idi, idj] = h**2/12
                        
    return Amatrix
                     


def Stencil_L(i):
    
    n = 2**i+1
       
    Lmatrix = eye((n-2)**2)
    
    for idi in range(0, (Lmatrix.shape[0])):
        
 
            
        
        negh_set = negh(n,  idi)
        

                
        for idj in negh_set:
                   

                        if idi - idj == 0:
                            
                            Lmatrix[idi, idj]= 4
                            
                           
                            
                        elif abs(idi - idj) == 1:
                            
                            Lmatrix[idi, idj] = -1
                        
                        elif abs(idi - idj) == int(n-2):
                            
                            Lmatrix[idi, idj] = -1

                        
    return Lmatrix




def FDM_G1(i):
    
    n = 2**i+1
    
    h =1/float(n-1)
       
    G1matrix = eye((n-2)**2)
    
    for idi in range(0, (G1matrix.shape[0])):
        
 
            
        
        negh_set = negh(n,  idi)
        

                
        for idj in negh_set:
                   

                    if idi-idj == 0:
                            
                            G1matrix[idi,idj] = 2*h
                            
                    elif abs(idi -idj) == int(n-2):
                            
                            G1matrix[idi,idj] = -h

    return G1matrix



def FDM_G2(i):
    
    n = 2**i+1
    
    h =1/float(n-1)
       
    G2matrix = eye((n-2)**2)
    
    for idi in range(0, (G2matrix.shape[0])):
        
 
            
        
        negh_set = negh(n,  idi)
        

                
        for idj in negh_set:
                   
    
                    if idi-idj == 0:
                            
                            G2matrix[idi,idj] = 2*h
                            
                    elif abs(idi -idj) == 1:
                            
                            G2matrix[idi,idj] = -h

    return G2matrix
    

def S_Preconditioner(i, alpha):
    
    import math
    

    
    
    Ama = Stencil_A(i)
    print Ama.shape, type(Ama)
    Lma = Stencil_L(i)
    print Lma.shape
    G1ma = FDM_G1(i)
    print G1ma.shape
    G2ma = FDM_G2(i)
    print G2ma.shape
#    ZeroMatrix = csr_matrix(((n-2)**2, (n-2)**2))
#    print ZeroMatrix.shape
#    

    
    
    SymBigMat = bmat([[Ama, None, None,  math.sqrt(alpha)*Lma],\
                       [None, Lma, None, G1ma.T],\
                       [None, None, Lma,  G2ma.T],\
                       [ math.sqrt(alpha)*Lma, G1ma, G2ma, None]])
    return SymBigMat
    



    