#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 14:12:46 2019

@author: shilu
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:46:38 2019

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
from scipy.sparse.linalg import spsolve,inv
from scipy.sparse.linalg import LinearOperator
from scipy import *
from Uniform_LinearOperator import Setup_LinearOperators, S_Preconditioner
from UniformTest import Setup_Matrices
from MG_Precon import multigrid_matvec, pre_multigrid
#from MG_uniform import rhs




maxIterations = 1000
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
    return 0

def QMR(A, left_precond, right_precond):
    #b = np.array([i for i in range(A.shape[0])])
    #b= np.array([[(i) for i in range(pre_fem_ini.shape[0])]]).T
    b=np.ones(A.shape[0])
    #rhs_list = np.reshape(rhs, (4*31*31,))
    #b = rhs_list
    #if (y == 0):
        
        #u = scipy.sparse.linalg.qmr(A,b,  tol=1e-6, maxiter=maxIterations,callback=Residual)
        
        #return u
        #
    #else:
    u = scipy.sparse.linalg.qmr(A,b,  tol=1e-6, maxiter=maxIterations, M1=left_precond, M2 = right_precond, callback=Residual)
        
    return u

# x=scipy.sparse.linalg.cg(A40520,np.ones(40520),maxiter=maxIterations)
        
    
def preCondition (A):
    from scipy.sparse.linalg import spilu
    Ainv=spilu(A,drop_tol= 1e-9)
    # Using invere.solve to define the mv of preconditioner  
    #hermitian_inv = spilu(A.T, drop_tol= 1e-9)
    
    preconditioner=LinearOperator(Ainv.shape, matvec=Ainv.solve, rmatvec=Ainv.solve)
    return preconditioner


def preCondition1(A):
    from scipy.sparse import csc_matrix
    from scipy.sparse.linalg import  spilu, splu#drop_tol= 1e-9)
    
    lu = splu(A)
    Pr = csc_matrix(A.shape)
    
    Pr[lu.perm_r, np.arange(A.shape[0])] = 1
    
    Pc = csc_matrix(A.shape)
    Pc[np.arange(A.shape[0]),lu.perm_c] = 1
    Left_pre = Pr.T* lu.L*lu.U*Pc.T
    hermitian =  (Left_pre).T
    
    Left_pre_inv = splu(Left_pre)# drop_tol= 1e-12)
    hermitian_inv = splu(hermitian)# drop_tol= 1e-12)
    
    Left_pre_inv=LinearOperator(Left_pre.shape, matvec=Left_pre_inv.solve , rmatvec=hermitian_inv.solve)
    return Left_pre_inv



def Left_preCondition(P):
    M_x = lambda x: spsolve(P, x)
    rM_x = lambda x: spsolve(P.T, x)
    M =  LinearOperator(P.shape, matvec= M_x, rmatvec= rM_x)
    return M
    
    
def conjGradError(A):
    b=np.ones(A.shape[0])
    return np.linalg.norm(b-A*QMR(A,preCondition(A),0)[0])/np.linalg.norm(b)




#fem = Setup_LinearOperators(i,alpha)
#pre_fem = Setup_Matrices(alpha, i, 'FEM')



#u_MG = QMR(fem, rhs_int_list,  pre_multigrid(n), preCondition(np.identity(4*n**2)))
#u_precond = QMR(fem, b, preCondition1(pre_fem), preCondition(np.identity(4*n**2)))
#u_precond = np.reshape(u_precond[0], (4,n,n))
#u_direct = QMR(fem, rhs_int_list, preCondition(np.identity(4*n**2)), preCondition(np.identity(4*n**2)))
#u_direct = np.reshape(u_direct[0], (4,n,n))

#print 'START'
#fem = Setup_LinearOperators(5,alpha)
#print 'FINISH'
#
#print 'START'
#pre_fem = S_Preconditioner(6,1e-09)
#print 'FINISH'


#b= np.array([[(i) for i in range(pre_fem_ini.shape[0])]]).T
#b = np.array([i for i in range(pre_fem_ini.shape[0])])
#b=np.ones(pre_fem_ini.shape[0])
#x0 = spsolve(pre_fem_ini, b)

fem, pre_fem = Setup_Matrices(1e-8, 5, 'FEM')

#pre_fem = Setup_Matrices(alpha, 5, 'FDM')


QMR(fem,preCondition1(pre_fem), preCondition1(np.identity(fem.shape[0])))

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

