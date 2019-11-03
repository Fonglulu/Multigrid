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
from UniformTest import Setup_Matrices, Setup_Rhs
from MG_Precon import multigrid_matvec, pre_multigrid
#from MG_uniform import rhs

#print(os.getcwd())

#Smatrix6 = np.load('Smatrix6.npy')
#Smatrix4_C = sparse.csr_matrix(Smatrix4_C).tocsc()
#femSmatrix6 = np.load('femSmatrix6.npy')
#femSmatrix4_C = sparse.csr_matrix(femSmatrix4_C).tocsc()

#os.chdir('/Users/shilu')
#print(os.getcwd())
#B162 = mmio.mmread("smdsurv162.mtx").tocsc()
#A40520 = mmio.mmread('smdsurv40520.mtx').tocsc()
#Atilde = mmio.mmread("modsurv40520.mtx").tocsc()
#Smatrix5 = np.load('Smatrix5.npy')
#Smatrix5 = sparse.csr_matrix(Smatrix5).tocsc()
#femSmatrix5 = np.load('femSmatrix5.npy')
#femSmatrix5 = sparse.csr_matrix(femSmatrix5).tocsc()


#print(type(femSmatrix5),type(Atilde))




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


def QMR(A,b, left_precond,right_precond):
    #b = np.array([i for i in range(A.shape[0])])
    #b= np.array([[(i) for i in range(pre_fem_ini.shape[0])]]).T
    #b=np.ones(A.shape[0])
    #rhs_list = np.reshape(rhs, (4*31*31,))
    #b = rhs_list
    #if (y == 0):
        
        #u = scipy.sparse.linalg.qmr(A,b,  tol=1e-6, maxiter=maxIterations,callback=Residual)
        
        #return u
        #
    #else:
    u = scipy.sparse.linalg.qmr(A,b,  tol=1e-06,  maxiter=maxIterations, M1=left_precond, M2=right_precond ,callback=Residual)
        
    return u

# x=scipy.sparse.linalg.cg(A40520,np.ones(40520),maxiter=maxIterations)
        
    
def preCondition (A):
    from scipy.sparse.linalg import spilu
    Ainv=spilu(A,drop_tol= 1e-6)
    # Using invere.solve to define the mv of preconditioner  
    #hermitian_inv = spilu(A.T, drop_tol= 1e-9)
    
    preconditioner=LinearOperator(Ainv.shape, matvec=Ainv.solve, rmatvec=Ainv.solve)
    return preconditioner


def preCondition1(A):
    from scipy.sparse import csc_matrix,lil_matrix
    from scipy.sparse.linalg import  spilu, splu#drop_tol= 1e-9)
    
    lu = splu(A)
    Pr = lil_matrix(A.shape)
    
    Pr[lu.perm_r, np.arange(A.shape[0])] = 1
    
    Pc = lil_matrix(A.shape)
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

#
#alpha = 1
#i = 5
#fem = Setup_LinearOperators(i,alpha)
#pre_fem = Setup_Matrices(alpha, i, 'FDM')[1]
#n = 2**i-1
#b=np.ones(pre_fem.shape[0])

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

#fem = Setup_Matrices(1, 5, 'FEM')
#pre_fem = Setup_Matrices(alpha, 5, 'FDM')




def plot_residual(alpha,i, grid_type):
    
    global residualVector
    from numpy import linalg as LA
    import numpy as np
    from numpy.linalg import inv
    from PlotRountines import smoother_plot
    
    
    
    residualVector = np.array([])
    
    if grid_type == 'odd':
        
        n = 2**i-1
    else:
        
        n= 2**i-2
    #fem = Setup_LinearOperators(i,alpha)
    fem = Setup_Matrices(alpha, i, 'FEM')[0]
    
    pre_fem = Setup_Matrices(alpha, i, 'FDM')[0]
    
    #pre_fem_ini = Setup_Matrices(1,i, 'FEM')[0]

    

    
    b=Setup_Rhs(alpha,i)[0]
    
    #b_ini = Setup_Rhs(1,i)
    
    
    #ini = spsolve(pre_fem_ini, b_ini)
#    ini_plot = np.reshape(ini[0:n**2], (n,n))
#    approx_plot(ini_plot)
    
    
    
    u_direct = spsolve(fem, b)
    
    
    
    u_direct = np.reshape(u_direct[0:n**2], (n,n))
    
    approx_plot(u_direct) 

    
    #print LA.cond(np.matmul(np.matrix(inv(pre_fem)), np.matrix(fem))), fem.shape[0], b.shape[0]

    #pre_multigrid(n)
    #QMR(fem, pre_multigrid(n), preCondition(np.identity(4*n**2)))
    u = QMR(fem,b,  preCondition1(pre_fem), preCondition1(np.identity(4*n**2)))

    #QMR(femSmatrix2, preCondition1(Smatrix2),preCondition(np.identity(femSmatrix2.shape[0])),1)
    u = np.reshape(u[0][0:n**2], (n,n))
         
    approx_plot(u)                                         
                                                    
                                                    
    #QMR(femSmatrix5_A, Left_preCondition(Smatrix5_A),preCondition(np.identity(femSmatrix5_A.shape[0])),1)
    #QMR(femSmatrix_D_sparse, Left_preCondition(Smatrix5_D_sparse), preCondition(np.identity(Smatrix5_A.shape[0])),1)
    #QMR(femSmatrix5, Left_preCondition(Smatrix5), preCondition(np.identity(femSmatrix5.shape[0])),1)
    #QMR(femSmatrix3_C, Left_preCondition(Smatrix3_C),preCondition(np.identity(femSmatrix3_C.shape[0])),1)
    i = 0
#                                                        
#    print len(residualVector)
#    while (len(residualVector) < maxIterations):
#        residualVector = np.append(residualVector,0)
#        i += 1
#                                                            
#        print i
    print residualVector, len(residualVector)
    plt.figure(1)
    plt.title('QMR Residual Plot for alpha =1e-9')
    plt.xlabel("Iteration")
    plt.ylabel("Residual")
    plt.semilogy(range(len(residualVector)), residualVector)# linewidth=0.3,c='b')
    plt.savefig('CondResPlotAtole7.png', format='png', dpi=500)
    plt.show()
    plt.clf()




def approx_plot(u):
    
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    #from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np

    
    # Get number of inteior nodes
    [xdim, ydim] = u.shape
    
    fig = plt.figure(figsize=(8,5))
    ax = fig.gca(projection='3d')
    
    entire_u = np.zeros((xdim+2, ydim+2))
    
    entire_u[1:-1,1:-1] = u
    
    
    # Get the spacing, for example, 5 nodes on side yNlields the spacing h =1/4
    h=1/float(xdim+1)
    
    x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    


    ax.plot_surface(x1, y1, entire_u ,cmap='viridis',linewidth=0)

    
    # Make the ticks looks pretty
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    plt.show()    
    
