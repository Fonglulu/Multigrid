#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 16:30:11 2020

@author: shilu
"""

import numpy as np
from numpy import linalg as LA

A = np.array([[2,2,3],[0,2,3],[2,0,2]])

x = np.array([[0.74958822, 0.45763206, 0.47821585]])

v= np.array([[3, 3, 3]])

B = A-np.matmul(x.T,v)
def Relaxation_Jacobi(A, b, num_itr):
    
    size = np.shape(A)[0]
    
#    x_new = np.zeros(size)
    x_new = np.array([1000,180000,60000])
    
    diag = np.diag(A)
    
    d = np.array([1/float(i) for i in diag])
    
    invD = np.diag(d)
    print invD
    
#    print (LA.eig(np.matmul(invD, A)))
    print (LA.eig(A))
    
    for i in range(num_itr):
        
#        x_new = np.zeros(size)
        
        rhs = b - np.matmul(A,x_new)
#        print x_new
        x_new = x_new + 0.6*rhs
#        print  x_new
    return x_new
    
    