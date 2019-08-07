#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 16:31:17 2018

@author: shilu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 09:43:17 2018

@author: FENG Shi Lu
"""

import numpy as np
from numpy import pi, sin, cos, exp, inf
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from pylab import title, xlabel, ylabel, clf, plot,show, legend
np.set_printoptions(precision=2)




###############################################################################




def exp_soln(x,y):
    """Define the function exp(x_0)exp(x_1)"""
    return exp(x)*exp(y)
    
def exp_rhs(x,y):
    """Define the function -2*exp(x_0)exp(x_1)"""
    return -2.0*exp(x)*exp(y)

def sin_soln(x,y):
    """Define the function sin(pi x_0)sin(pi x_1)"""
    return sin(pi*x)*sin(pi*y)
    
def sin_rhs(x,y):
    """Define the function 2pi^2 sin(pi x_0)sin(pi x_1)"""
    return 2.0*pi*pi*sin(pi*x)*sin(pi*y)



###############################################################################
    








def Injection(uf):
    
    """ 
    Restrict find grid to coarse grid by injection
    
    Input: current approximation on fine grid uf
    
    Output: Restricted approximation on coarse grid 
    
    """
    
    # Get the current size of approximation
    [xdim, ydim] = uf.shape
    
    # Restrict on coarse grid
    return uf[0:xdim:2,0:ydim:2]
    

def Interpolation(uc):
    
    """ 
    Interpolate coarse grid to fine grid
    
    Input: current approximation on coarse grid uc

    Output: Interpolated approximation on find grid    
    """
    #print uc, 'BEfore'
    # Get the current size of approximation on coarse grid
    [xdim, ydim] = uc.shape
    
    # Initialise a next fine grid
    xnodes = 2*xdim-1
    ynodes = 2*ydim-1
    grid = np.zeros((xnodes,ynodes))
    
    
    # For even ordered i and j
    for i in range(xdim):
        for j in range (ydim):
            grid[2*i,2*j]=uc[i,j]
    

    # For even ordered j     
    for i in range(0, ynodes, 2):
        for j in range(1, xnodes-1, 2):
            grid[i,j]=0.5*(grid[i,j-1]+grid[i,j+1])

        
    # For even ordered i   
    for i in range(1, xnodes-1, 2):
        for j in range (0, ynodes, 2):
            grid[i,j]=0.5*(grid[i-1,j]+grid[i+1,j])
    
    # For odd ordered i and j
    for i in range (1, xnodes-1, 2):
        for j in range (1, ynodes-1, 2):
            grid[i,j]=0.25*(grid[i-1,j]+grid[i+1,j]+grid[i,j-1]+grid[i,j+1])
    
    return grid



###############################################################################
def Jacobi(u,rhs,n):
    """ 
    Jacobi iterative method
    Input: u: current approximation 
           rhs: rhs load vector
           n: grid size
           
    Output: new approximation u
    """       
    
    # Get the current size of RHS function 
    [xdim,ydim] = rhs.shape
    
    # Get the spacing
    h=1/float(n-1)
    
    # store the current approximation u as v
    v=u
    
    # use v to find the new approximation
    for i in range(1,xdim-1):
        for j in range(1,ydim-1):
            u[i,j]= (v[i-1,j]+v[i+1,j]+v[i,j-1]+v[i,j+1]+(h**2)*rhs[i,j])/4
    return u



def GaussSeidel(u,rhs,n):
    """
    Gauss-Seidel iterative method
    Input: u: current approximation
           rhs: rhs load vector
           n: grid size
           
    Output: new approximation u
    """
    #print u, 'ubefore'
    # Get the current size of RHS function
    [xdim,ydim] = rhs.shape
    
    # Get the spacing
   # print n, 'n'
    h=1/float(n-1)
    
    # Immeidately use the new approximation for next iteration
    for i in range(1,xdim-1):
        for j in range(1,ydim-1):
            u[i,j]= (u[i-1,j]+u[i+1,j]+u[i,j-1]+u[i,j+1]+(h**2)*rhs[i,j])/4
    
    #print u, 'uafter'
    return u
###############################################################################    



def residue(rhs, u, n):
    
    """ 
    Compute the residual r = f - Av
    Input: rhs: RHS function
           u: current approximation
           n: grid size
    
    Output: the residual r
    """
    
    
    # Get the current size of RHS function
    [xdim,ydim] = rhs.shape
    
    # Initialise the residual
    r=np.zeros((xdim,ydim))
    h=1/float(n-1)

    # Find each component of residual
    for i in range(1, xdim-1):
        for j in range(1, xdim-1):
            r[i,j] = rhs[i,j]+ (u[i-1,j]+u[i+1,j]+u[i,j-1]+u[i,j+1]-4*u[i,j])/(h**2)
    return r
    
###############################################################################
#    
# Input: RHS function of Poisson equation
#      : u initial guess
#      : rhs function
def VCycle(rhs, u, n, s1, s2, s3):
    """ Conduct one V-cycle for Poisson equation
        Input: rhs: RHS function of Poisson equation
               u:  initial guess
               n: initial mesh size n
               s1: unward sweeps  
               s2: coareset sweeps
               s3: Upward sweeps
               
        Output: Improved inital guess after a cycle
    """
    
    # Set up the finest grid
    
    # Mesh spacing on finest grid
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
    
    # Initialise a list records the approximation on different level
    u_list = []
    
    # x directio of the finest grid
    xdim.append(rhs.shape[0])
    
    # y direction of the finest grid
    ydim.append(rhs.shape[1])
    
    # h avlue of the finest grid
    hvalue.append(h)
    
    # rhs on finest grid
    rhs_list.append(rhs)
    
    # u on finest grid
    u_list.append(u)

    
    
    
    # Set up the coarse grids
    
    
    # the numer nodes on coarse grid is rougly a half from upper fine grid
    while( (xdim[levelnumber]-1) % 2 == 0 and (ydim[levelnumber]-1) % 2 ==0\
          
           and xdim[levelnumber]-1 >2 and ydim[levelnumber]-1 >2 ):
        
        
        levelnumber = levelnumber+1
        
        xdim.append((xdim[levelnumber - 1] -1) //2 +1)
        
        ydim.append((ydim[levelnumber - 1] -1) //2 +1)
        
        hvalue.append(2*hvalue[levelnumber-1])
        
        # Append rhs for following  meshs as temporarily putting zero vector
        rhs_list.append(np.zeros((xdim[levelnumber],ydim[levelnumber])))
        
        # Append u for following  meshs as temporarily putting zero vector
        u_list.append(np.zeros((xdim[levelnumber],ydim[levelnumber])))
        
        
    totallevel = levelnumber
    
    #print totallevel
        
        


    # Downward 
    
    
    for levelnumber in range(totallevel):

        # Get the initialised u from each level
        ulevel = u_list[levelnumber]
        
        
        
        
        
        #print(u_list[0])
 ##############################################################################       
#        fig = plt.figure(figsize=(8,5))
#        ax = fig.gca(projection='3d')
#        
#        x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
#    
#    
#    #uexact = sin_soln(x1, y1)
#    # Plot the surface.
#        surf = ax.plot_surface(x1, y1, u_list[0], cmap=cm.coolwarm,
#                           linewidth=0, antialiased=False)
#    
#    # Customize the z axis.
#    #ax.set_zlim(-1.01, 1.01)
#        ax.zaxis.set_major_locator(LinearLocator(10))
#        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#    
#    # Add a color bar which maps values to colors.
#        fig.colorbar(surf, shrink=0.5, aspect=5)
#        plt.show()
##############################################################################        
        
        
        
        # Get the initialised rhs from each level
        rhslevel = rhs_list[levelnumber]
        
        #  Apply s1 times smooter on each level  to solve Au =f (r)
        for sweeps in range(s1):
            
            ulevel = GaussSeidel(ulevel, rhslevel, xdim[levelnumber])

        u_list[levelnumber]=ulevel
        
 ##############################################################################       
#        fig = plt.figure(figsize=(8,5))
#        ax = fig.gca(projection='3d')
#        
#        x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
#    
#    
#    #uexact = sin_soln(x1, y1)
#    # Plot the surface.
#        surf = ax.plot_surface(x1, y1, u_list[0], cmap=cm.coolwarm,
#                           linewidth=0, antialiased=False)
#    
#    # Customize the z axis.
#    #ax.set_zlim(-1.01, 1.01)
#        ax.zaxis.set_major_locator(LinearLocator(10))
#        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#    
#    # Add a color bar which maps values to colors.
#        fig.colorbar(surf, shrink=0.5, aspect=5)
#        plt.show()
############################################################################## 
        
        
        
        # Get the residual after relax
        r = residue (rhslevel, ulevel, xdim[levelnumber])
        
        # Replace the old rhs vector by new residual
        rhs_list[levelnumber+1] = Injection(r)






        
    
    # Solve for the error on coarest grid
    
    # Find the u on coarest grid
    ulevel = u_list[totallevel]
    
    # Find the rhs on coarest grid 
    rhslevel = rhs_list[totallevel]
    
    # Apply s2 times smoother on coarest grid to solve Au =f (r)
    for sweeps in range(s2):
        
        ulevel = GaussSeidel(ulevel, rhslevel, xdim[totallevel])
    
    
    # Replace the old approximation on coarest grid by new one
    u_list[totallevel] = ulevel
    
    

    # Upward
    
    for levelnumber in range(totallevel-1, -1, -1):
        
        # For each coarse grid
        udown = u_list[levelnumber+1]
        
        # Find the upper fine grid from initilsed u list
        ulevel = u_list[levelnumber]
        
        # Update the new approximation u on the fine grid
        ulevel = ulevel + Interpolation(udown)
        
        # Get the rhs vector on each level
        rhslevel = rhs_list[levelnumber]
        
        # Apply s3 times smoother on fine grid to solve Au =f (r)
        for sweeps in range(s3):
            
            ulevel = GaussSeidel( ulevel, rhslevel, xdim[levelnumber])
            
        r = residue( rhslevel, ulevel, xdim[levelnumber])
        
        # Update the approximation on each level
        u_list[levelnumber] = ulevel
        
      
        
    # Update the improved initial guess    
    u=u_list[0]
    
    return u
    
###############################################################################

def MGsolverUP(cyclenumber, rhs, exact, levelnumber):
    
    # Set up the grid size
    #print rhs, 'RRHHSS'
    n= 2**(levelnumber+1)+1
    #print n, 'n'
    h=1/float(n-1)
    
    #print h

    
    # Set up the mesh 
#    x = np.linspace(0,1+h,h)
#    y = np.linspace(0,1+h,h)
#
#    x1, y1 = np.meshgrid(x,y)
#    print x1
    x1, y1 = np.meshgrid(np.arange(0, 1+h, h), np.arange(0, 1+h, h))
    #print x1
    
    # Set the initial guess
    #u=np.zeros((n,n))
    u= exact
    #print u.shape[0],' ushape'
    
    # Set the boundary (also the exact solution in this case)
    #uexact = zero(x1,y1)
    #print uexact
    
    # Set the rhs function
    #rhs = exp_rhs(x1,y1)

    #Set Boundary 
#    u[0] = exact(x1, y1)[0]
#    
#    u[-1] = exact(x1, y1)[-1]
#    
#    u[:,0]= exact(x1, y1)[:,0]
#    
#    u[:,-1] = exact(x1, y1)[:,-1]
    
    #print u, 'u'
    
    # Set the number of relax
    s1=2
    s2=2
    s3=2
    
    # Initialise a list to record l2 norm of resudual 
    rnorm=[np.linalg.norm(residue(rhs, u, n)) * h]
    
    # Initialise a list to record l2 norm of error
    #enorm=[np.linalg.norm(u-exact(x1,y1))*h]
    
    
    # Start V-cycle
    for cycle in range(1, cyclenumber+1):
        
        
        u = VCycle(rhs, u, n, s1, s2, s3)
        

        rnorm.append(np.linalg.norm(residue(rhs, u, n))*h)
        
        
        #enorm.append(np.linalg.norm(u-exact(x1,y1))*h)
        
     #Plot the semi-log for resiudal and erro
#    xline = np.arange(cyclenumber+1)
#    plt.figure(figsize=(4,5))
#    plt.semilogy(xline, enorm, 'bo', xline, enorm, 'k',label='sdad')
#    title('Convergence wrt Error (GaussSeidel)')
#    xlabel('Number of cycles')
#    ylabel('Error under l2 norm')
#    plt.show()
    #print rhs,'AFTERRHS'
    return u
