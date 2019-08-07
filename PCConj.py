#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 10:45:28 2019

@author: shilu
"""

import timeit
import inspect
import numpy as np
import scipy 
from scipy import sparse
import matplotlib.pyplot as plt
from scipy.io import mmio
import os

print(os.getcwd())

os.chdir('/Users/shilu')
print(os.getcwd())
#B162 = mmio.mmread("smdsurv162.mtx").tocsc()
A40520 = mmio.mmread('smdsurv40520.mtx').tocsc()
Atilde = mmio.mmread("modsurv40520.mtx").tocsc()


maxIterations = 300
iterSpace = np.linspace(0,maxIterations,num=maxIterations)
residualVector = np.array([])
sparseResidual= []

def sparseDirect (A):
    b=np.ones(A.shape[0])
    return scipy.sparse.linalg.spsolve(A,b)

def sparseError(A):
    b=np.ones(A.shape[0])
    return np.linalg.norm(A*sparseDirect(A)-b)/np.linalg.norm(b)

def Residual (xk):
    frame = inspect.currentframe().f_back
    global residualVector
    residualVector = np.append(residualVector,frame.f_locals['resid'])
    return -1


def conjGrad(A,precond,y):
    b=np.ones(A.shape[0])
    #if (y == 0):
        #return scipy.sparse.linalg.cg(A,b,tol=1e-3,maxiter=maxIterations,callback=Residual)
    #else:
    return scipy.sparse.linalg.cg(A,b,tol=1e-3, maxiter=maxIterations,callback=Residual,M=precond)

# x=scipy.sparse.linalg.cg(A40520,np.ones(40520),maxiter=maxIterations)
def preCondition (A):
    from scipy.sparse.linalg import spilu
    Ainv=spilu(A,drop_tol= 1e-9)
    
    # Using invere.solve to define the mv of preconditioner  
    preconditioner=scipy.sparse.linalg.LinearOperator(Ainv.shape,matvec=Ainv.solve)
    return preconditioner

def conjGradError(A):
    b=np.ones(A.shape[0])
    return np.linalg.norm(b-A*conjGrad(A,preCondition(A),0)[0])/np.linalg.norm(b)


#conjGrad(Atilde,preCondition(np.identity(Atilde.shape[0])),1)
conjGrad(Atilde,preCondition(Atilde),1)

i =0 

print len(residualVector)
while (len(residualVector) < maxIterations):
    i +=1
    print i
    residualVector = np.append(residualVector,0)
    
    
   

plt.figure(1)
plt.title('Conjugate Gradient Residual Plot for A (conditioned 1e-9)')
plt.xlabel("Iteration")
plt.ylabel("Residual")
plt.semilogy(iterSpace, residualVector)#linewidth=0.3,c='b',s =1)
plt.savefig('CondResPlotAtole7.png', format='png', dpi=500)
plt.show()
plt.clf()

#def plot_residue(Atilde, matrix):
#
#    conjGrad(Atilde,preCondition(matrix),1)
#    
#    i= 0
#
#    while (len(residualVector) < 100):
#        print len(residualVector)
#        
#        
#        
#        residual = np.append(residualVector,0)
#        
#        i+=1
#        print(i)
#        
#    plt.figure(1)
#    plt.title('Conjugate Gradient Residual Plot for A (conditioned 1e-9)')
#    plt.xlabel("Iteration")
#    plt.ylabel("Residual")
#    plt.scatter(iterSpace, residual, linewidth=0.3,c='b')
#    plt.savefig('CondResPlotAtole7.png', format='png', dpi=500)
#    plt.show()
#    plt.clf()