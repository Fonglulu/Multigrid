#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 19:57:21 2018

@author: FENG Shi Lu

This module stores the smoothers
"""


import numpy as np
from numpy import  pi, sin, cos, exp, inf
from scipy.sparse.linalg import spsolve
from scipy import eye, zeros, linalg
from numpy import linalg as LA

from JacobiDown import JacobiDown
#from JacobiUp import JacobiUp
#from MGsolver import MGsolve
#from MGsolverUp import MGsolverUP

from grid.Grid import Grid
from grid.NodeTable import not_slave_node
from BuildSquare import build_square_grid
from BuildEquation import build_equation_linear_2D, set_polynomial_linear_2D,\
Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY, NIn_triangle, build_matrix_fem_2D
from grid.function.FunctionStore import zero, linear

from operator import itemgetter
np.set_printoptions(precision=4)

from copy import copy


def sin_soln(x,y):
    """Define the function sin(pi x_0)sin(pi x_1)"""
    return sin(pi*x)*sin(pi*y)

def cos_sin(x, y):
    """Define the function pi cos(pi x_1)sin(pi y_1)"""
    
    return pi*cos(pi*x)*sin(pi*y)

def sin_cos(x, y):
    """Define the function pi sin(pi x_1)cos(pi y_1)"""
    
    return pi*sin(pi*x)*cos(pi*y)

def exp_soln(x,y):
    """Define the function exp(x_0)exp(x_1)"""
    return exp(x)*exp(y)
    

def exp_w(x,y):
    """Define the function -2*exp(x_0)exp(x_1)"""
    return -1.0*exp(x)*exp(y)


def Sin4(x,y):
    
    return sin(4*pi*x)*sin(4*pi*y)

def Xsin4(x,y):
    
    return 4*pi*cos(4*pi*x)*sin(4*pi*y)

def Ysin4(x,y):
    
    return 4*pi*sin(4*pi*x)*cos(4*pi*y)

def XYsin4(x,y):
    
    return - (-1*32*pi*pi*sin(4*pi*x)*sin(4*pi*y))

def plain(x,y):
    
    return x-x+1

def Zero(x,y):
    
    return x-x+0.0

def Linear(x,y):
    
    return x+y

def Xlinear(x,y):
    
    return x-x+1.0
    

def Ylinear(x,y):
    
    return y-y+1.0


def L3(x,y):
    
    return x**3 + y**3

def x_3(x,y):
    
    return 3*x**2


def y_3(x,y):
    
    return 3*y**2 

def XYl_3(x,y):
    
     return -6*0.01*(x+y)
    



def Injection(uf):
    
    """ 
    Restrict find grid to coarse grid by injection
    
    Input: current approximation on fine grid uf
    
    Output: Restricted approximation on coarse grid 
    
    """
    #print uf[1], 'injectbefore'
    # Get the current size of approximation
    [depth, xdim, ydim] = uf.shape
    
    # Restrict on coarse grid
    
    return uf[:, 0:xdim:2, 0:ydim:2]



def Interpolation(uc):
    
    """ 
    Interpolate coarse grid to fine grid
    
    Input: current approximation on coarse grid uc

    Output: Interpolated approximation on find grid    
    """
    [depth, xdim, ydim] = uc.shape
    #print depth, xdim, ydim
    
    # Initialise a next fine grid
    xnodes = 2*xdim-1
    ynodes = 2*ydim-1
    grid = np.zeros((depth, xnodes,ynodes))
    
    
    # For even ordered i and j
    for k in range(depth):
        for i in range(xdim):
            for j in range (ydim):
                grid[k, 2*i, 2*j]=uc[k, i,j]
    

    # For even ordered j  
    for k in range(depth):
        for i in range(0, ynodes, 2):
            for j in range(1, xnodes-1, 2):
                grid[k,i,j]=0.5*(grid[k,i,j-1]+grid[k,i,j+1])

        
    # For even ordered i   
    for k in range(depth):
        for i in range(1, xnodes-1, 2):
            for j in range (0, ynodes, 2):
                grid[k,i,j]=0.5*(grid[k,i-1,j]+grid[k,i+1,j])
    
    # For odd ordered i and j
    for k in range(depth):
        for i in range (1, xnodes-1, 2):
            for j in range (1, ynodes-1, 2):
                grid[k,i,j]=0.25*(grid[k,i-1,j]+grid[k,i+1,j]+grid[k,i,j-1]+grid[k,i,j+1])#    

            
            
    return grid





def residue(totalrhs,  u, alpha):
    
    #print totalrhs[0], 'RRS'\
    # Get the current size of RHS function
    [depth, xdim,ydim] = totalrhs.shape
    #print 'rhs', rhs
    
    
    # Initialise the residual
    r=np.zeros((depth, xdim,ydim))
    
    h = 1/ float(xdim -1)
    
    for i in range(1, xdim-1):
        for j in range(1, xdim-1):
            
            # rhs[0] = f1 +G1g1 + G2g2
            r[0,i,j] = totalrhs[0, i,j] +(u[0, i-1,j]+u[0, i+1,j]+u[0, i,j-1]+u[0, i,j+1]-4*u[0, i,j])
            

            r[1,i,j] = totalrhs[1, i,j] + alpha*(u[1, i-1,j]+u[1, i+1,j]+u[1, i,j-1]+u[1, i,j+1]-4*u[1, i,j])
            
            r[2,i,j] = totalrhs[2, i,j]+ alpha* (u[2, i-1,j]+u[2, i+1,j]+u[2, i,j-1]+u[2, i,j+1]-4*u[2,i,j])
            
            r[3,i,j] = totalrhs[3, i,j] + (u[3, i-1,j]+u[3, i+1,j]+u[3, i,j-1]+u[3, i,j+1]-4*u[3,i,j])
            
    #print r[1], 'r'
    
    
    
    
    #print r[0], 'Residual'
    return r









def VCycle(total, rhs, u, n, s1, s2, s3, dataX, dataY,  data, alpha):
    
    # Set up the finest grid
    
    # Mesh spacing on finest grid, for debugging
    h = 1/float(n-1)
    
    # Levelnumber of finiest grid
    levelnumber = 0
    
        
    # Initialise a list records the x direction of grid
    xdim = []
    
    # Initialise a list records the y direction of grid
    ydim = []
    
    # Initialise a list records the spacing values
    hvalue = []
    
    # Initialise a list records the rhs on different level
    rhs_list = []
    
    rhs_total = []
    
    # Initialise a list records the approximation on different level
    u_list = []
    
    # x direction of the finest grid
    xdim.append(rhs.shape[1])
    #print xdim, 'lll'
    
    # y direction of the finest grid
    #print rhs.shape[1], rhs.shape[2], 'dadad'
    ydim.append(rhs.shape[2])
    #print ydim, 'yll'
    # h avlue of the finest grid
    hvalue.append(h)
    
    # rhs on finest grid
    rhs_list.append(rhs)
    
    rhs_total.append(total)
    
    # u on finest grid
    u_list.append(u)
    
    
    # Set up the coarse grids
    
    
    # the numer nodes on coarse grid is rougly a half from upper fine grid
    while( (xdim[levelnumber]-1) % 2 == 0 and (ydim[levelnumber]-1) % 2 ==0\
          
           and xdim[levelnumber]-1 >2 and ydim[levelnumber]-1 >2 ):
        
        #print 'i'
        levelnumber = levelnumber+1
        
        xdim.append((xdim[levelnumber - 1] -1) //2 +1 )
        #print xdim, 'ol'
        
        ydim.append((ydim[levelnumber - 1] -1) //2 +1)
        
        #hvalue.append(2*hvalue[levelnumber-1])
        
        # Append rhs for following  meshs as temporarily putting zero vector
        rhs_list.append(np.zeros((4, xdim[levelnumber],ydim[levelnumber])))
        
        # Append u for following  meshs as temporarily putting zero vector
        u_list.append(np.zeros((4, xdim[levelnumber],ydim[levelnumber])))
        
        rhs_total.append(np.zeros((4, xdim[levelnumber],ydim[levelnumber])))
        
        
    totallevel = levelnumber
    #print totallevel, 'dadas'
    
    for levelnumber in range(totallevel):

        # Get the initialised u from each level
        ulevel = u_list[levelnumber]
        
        # Get the initialised rhs from each level
        rhslevel = rhs_list[levelnumber]
        #print rhslevel, 'rhslevel'
        
        totalrhs = rhs_total[levelnumber]
            
        #  Apply s1 times smooter on each level  to solve Au =f (r)
        
    
        for sweeps in range(s1):
            
                
                ulevel = JacobiDown(ulevel, rhslevel, totalrhs, totallevel-levelnumber, dataX, dataY,\
                                data, alpha)

              
        u_list[levelnumber]=ulevel
        
        # Get the residual after relax
        #print totalrhs, 'totalrhs'
        r = residue (totalrhs,  ulevel, alpha)
        #print r, 'rhs'
        
        # Replace the old rhs vector by new residual
        rhs_list[levelnumber+1] = Injection(r)
        #print rhs_list
        
        
    # Solve for the error on coarest grid
    
    # Find the u on coarest grid
    ulevel = u_list[totallevel]
    
    #print ulevel[0], 'CORARSE'
    # Find the rhs on coarest grid 
    rhslevel = rhs_list[totallevel]
    
    
    
    totalrhs = rhs_total[totallevel]
    
    #print rhslevel[0], 'CORARSE'
    
    print 'Coarsest'
    # Apply s2 times smoother on coarest grid to solve Au =f (r)
    for sweeps in range(s2):
        
       
        ulevel = JacobiDown(ulevel, rhslevel, totalrhs,  totallevel-totallevel,  dataX, dataY, data, alpha)
    
    
    # Replace the old approximation on coarest grid by new one
    u_list[totallevel] = ulevel
    #print u_list, 'ULIST'
    
    print 'Upwards'
    
    # Upward
    
    for levelnumber in range(totallevel-1, -1, -1):
        
        # For each coarse grid
        udown = u_list[levelnumber+1]
        #print udown, 'udown'
        
        # Find the upper fine grid from initilsed u list
        ulevel = u_list[levelnumber]
        
        #print Interpolation(udown), 'inter'
        # Update the new approximation u on the fine grid
        ulevel = ulevel + Interpolation(udown)
        
        # Get the rhs vector on each level
        rhslevel = rhs_list[levelnumber]
        
        totalrhs = rhs_total[levelnumber]
        
        # Apply s3 times smoother on fine grid to solve Au =f (r)
        for sweeps in range(s3):
            
#            ulevel = JacobiUp( ulevel, rhslevel, totalrhs, totallevel-levelnumber, dataX, dataY,\
#                             data, alpha)
            ulevel = JacobiDown(ulevel, rhslevel, totalrhs,  totallevel-levelnumber, dataX, dataY, data, alpha)
        #r = residue( rhslevel, ulevel, alpha)
        
        # Update the approximation on each level
        u_list[levelnumber] = ulevel
        
      
        
    # Update the improved initial guess    
    u=u_list[0]
    
    return u
    

    
    
    

def Ultimate_MG(cyclenumber):
    
    from TPSFEM import Amatrix, h1, h2, h3, h4, Lmatrix, G1, G2, dvector
    
    # Set up the grid size， this includes boundary

    grid = Grid()
    
    i = 2

    n= 2**i+1
    
    # Find the spacing
    h=1/float(n-1)
    
    # Set the mesh grid, that is interior
    #x1, y1 = np.meshgrid(np.arange(h, 1, h), np.arange(h, 1, h))
    x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    #print x1, y1
    # Set the data space
    datax = np.linspace(0, 1.0,30)
    datay = np.linspace(0, 1.0,30)
  
    dataX, dataY = np.meshgrid(datax,datay)
    

    # Set the exact solution of c, g1, g2, w on every node
    cexact = L3
    
    g1exact =  x_3
    
    g2exact = y_3
    
    wexact = XYl_3
    
    # Set the shape of objective function
    Z = cexact(dataX, dataY)
    
    
    data = Z.flatten()
    

    # Set penalty term
    alpha = 0.01
    
    # Initialise grid for calculating matrices and boundaries
    grid = Grid()
    
    # Build square 
    build_square_grid(n, grid, zero)
    
    # Store matrices on grid
    build_matrix_fem_2D(grid, Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY,  dataX, dataY)
    
#    Nl=[]
#    for node in (not_slave_node(grid)):
#        
#        Nl.append([node,node.get_node_id()])
#        
#    Nl=sorted(Nl, key = itemgetter(1))
#    
#    Nl=[node[0] for node in Nl]
            
#    for node in Nl:
#        node.set_value(Nl.index(node))
        
    
    # Find the boundary for Ac + Lw = d - h1
    h1 = np.reshape(h1(grid,cexact, wexact), (n-2,n-2))
    #h1 = h1(grid,cexact, wexact)
    
    # Find the boundary for alphaLg1 -G1^Tw = -h2
    h2 = np.reshape(h2(grid, g1exact, wexact, alpha), (n-2,n-2))
    #h2 = h2(grid, g1exact, wexact, alpha)
    
    # Find the boundary for alphaLg2 -G2^Tw = -h3
    h3 = np.reshape(h3(grid, g2exact, wexact, alpha), (n-2,n-2))
    #h3 = h3(grid, g2exact, wexact, alpha)
    
    # Find the boundary for Lc -G1g1 -G2g2 = -h4
    h4 = np.reshape(h4(grid, cexact, g1exact, g2exact), (n-2,n-2))
    #h4 = h4(grid, cexact, g1exact, g2exact)
    
    # Find d vector
    #print dvector(grid, data)
    
    dvector = np.reshape(dvector(grid, data), (n-2,n-2))
    
    #dd = dvectorgrid, data
    
    #print dd
    

    

    
    # Set the initial guess for interior nodes values
    #u=np.zeros((4,n-2, n-2))
    u=200*np.ones((4,n,n))
    
    # Set RHS at intilisation
    #rhs = np.zeros((4, n-2, n-2))
    rhs = np.zeros((4,n,n))
    rhs[0][1:-1,1:-1] = -h4
    #print rhs[0][1:-1,1:-1]
    rhs[1][1:-1,1:-1] = -h2
    #print rhs[1], 'rhs[1]'
    rhs[2][1:-1,1:-1] = -h3
    rhs[3][1:-1,1:-1] = dvector-h1
    #print rhs[3][1:-1,1:-1]
    
    total = np.zeros((4,n,n))
    
    # Set the boundary (also the exact solution in this case)
    
    u[0,0,:] = cexact(x1, y1)[0]
    
    
    u[0, -1,:] = cexact(x1, y1)[-1]
    
    u[0, :, 0] = cexact(x1, y1)[:, 0]
    
    u[0,:,-1] = cexact(x1, y1)[:,-1]
    
    u[1,0,:] = g1exact(x1, y1)[0]
    
    u[1, -1,:] = g1exact(x1, y1)[-1]
    
    u[1, :, 0] = g1exact(x1, y1)[:,0]
    
    u[1,:,-1] = g1exact(x1, y1)[:,-1]
    
    u[2,0,:] = g2exact(x1, y1)[0]
    
    u[2, -1,:] = g2exact(x1, y1)[-1]
    
    u[2, :, 0] = g2exact(x1, y1)[:,0]
    
    u[2,:,-1] = g2exact(x1, y1)[:,-1]
    
    u[3,0,:] = wexact(x1, y1)[0]
    
    u[3, -1,:] = wexact(x1, y1)[-1]
    
    u[3, :, 0] = wexact(x1, y1)[:,0]
    
    u[3,:,-1] = wexact(x1, y1)[:,-1]
    
    # Set the number of relax
    s1=2
    s2=2
    s3=2
    
    
    #Initialise a list to record l2 norm of resudual 
    rnorm=[np.linalg.norm(residue(rhs, u, alpha)) * h]
    
    # Initialise a list to record l2 norm of error
    ecnorm = [np.linalg.norm(u[0]-cexact(x1, y1))*h]
#    e2norm = [np.linalg.norm(u[1]-g1exact(x1, y1))*h]
#    e3norm = [np.linalg.norm(u[2]-g2exact(x1, y1))*h]
#    e4norm = [np.linalg.norm(u[3]-cexact(x1, y1))*h]
    
    
    
    # Start V-cycle
    for cycle in range(1, cyclenumber+1):
        
        
        u = VCycle(total, rhs, u, n, s1, s2, s3, dataX, dataY, data, alpha)
    
        

        rnorm.append(np.linalg.norm(residue(rhs, u, n))*h)
        
        
        ecnorm.append(np.linalg.norm(u-cexact(x1,y1))*h)
        
    print rnorm, 'rnorm'   
    print ecnorm, 'ecnorm'
    return u
    
    

    # Set the RHS 
    
    
    
def Plot_approximation(cyclenumber):   
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np

    
    
    fig = plt.figure(figsize=(8,5))
    ax = fig.gca(projection='3d')
    
    # Make data.
    i = 2
    
    n= 2**i+1
    
    h=1/float(n-1)
    
    x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    
    #u = sin(pi*x1)*sin(pi*y1) + 5*sin(31*pi*x1)*sin(31*pi*y1)
    u = Ultimate_MG(cyclenumber)[0]
    print(u)
    
    #uexact = exp_soln(x1, y1)
     #Plot the surface.
#    surf = ax.plot_surface(x1, y1, u, cmap=cm.coolwarm,
#                           linewidth=0, antialiased=False)
#    
#    # Customize the z axis.
#    #ax.set_zlim(-1.01, 1.01)
#    ax.zaxis.set_major_locator(LinearLocator(10))
#    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#    
#    # Add a color bar which maps values to colors.
#    fig.colorbar(surf, shrink=0.5, aspect=5)
     
    ax.plot_surface(x1, y1, u,cmap='viridis',linewidth=0)
    
    # Set the z axis limits
    #ax.set_zlim(node_v.min(),node_v.max())
    
    # Make the ticks looks pretty
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    plt.show()    
    
    