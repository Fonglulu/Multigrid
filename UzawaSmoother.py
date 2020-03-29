#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:57:59 2019

@author: shilu
"""


import numpy as np
from scipy.sparse.linalg import spsolve




def Lstencil(u):
    
    """This routine builds the stencil for L matrix"""
    
    newu  = np.zeros((u.shape[0], u.shape[0]))
    
    for i in range(1, u.shape[0]-1):
        
        for j in range(1, u.shape[0]-1):
            
            #print u[i,j] ,u [i-1,j] ,u[i+1,j], u[i,j+1]
            newu[i,j] = (4* u[i,j] -  u[i-1,j] - u[i+1,j] - u[i,j-1] - u[i, j+1])
            
    return newu



def Astencil(u,h):
    
    """ This routine builds the stencil for A matrix"""
    
    newu  = np.zeros((u.shape[0], u.shape[0]))
    
    for i in range(1, u.shape[0]-1):
        
        for j in range(1, u.shape[0]-1):
            
            newu[i,j] = (6*u[i,j] +u[i,j+1]+u[i,j-1] + u[i-1,j] + u[i+1, j] +u[i+1,j+1] + u[i-1,j-1])*(h**2)/float(12)
            
    return newu


def G1stencil(u,h):
    
    """This routine builds the stencil for G1 matrix"""
    
    newu  = np.zeros((u.shape[0], u.shape[0]))
    
    for x in range(1, u.shape[0]-1):
        
        for y in range(1, u.shape[0]-1):
            
            # Swap x,y node-edge data structure grid
            #newu[x,y] = (-2*u[x,y-1]+2*u[x,y+1] - u[x+1, y] + u[x-1, y] +u[x+1, y+1] - u[x-1, y-1])*h/float(6)
            

            

             # correct G1 stencil
            newu[x,y] = (-2*u[x-1, y]+2*u[x+1,y] - u[x, y+1] + u[x, y-1] - u[x-1, y-1] +u[x+1, y+1])*h/float(6)
            
#            newu[x,y] = (-u[x, y] +u[x+1,y])*h
            
            # linda's version
            #newu[x,y] = (2*u[x,y-1]-2*u[x,y+1] - u[x+1, y] + u[x-1, y] -u[x+1, y+1] + u[x-1, y-1])*h/float(6)
            




            
    return newu


def G2stencil(u,h):
    
    """This rountine builds te stencil for G2 matrix"""
    
    newu  = np.zeros((u.shape[0], u.shape[0]))
    
    for x in range(1, u.shape[0]-1):
        
        for y in range(1, u.shape[0]-1):
            
            #Swap x,y node-edge data structure grid
            #newu[x,y] = (u[x,y-1] - u[x,y+1] +2 *u[x+1,y] -2*u[x-1,y] +u[x+1, y+1] -u[x-1,y-1])*h/float(6)
            
  
            # correct G2 stencil
            newu[x,y] = (u[x-1,y] - u[x+1,y] +2 * u[x,y+1] -2*u[x,y-1] - u[x-1, y-1] + u[x+1,y+1])*h/float(6)
#            newu[x,y] = (-u[x,y]+u[x,y+1])*h
            
            
            # linda's version
            #newu[x,y] = (u[x,y-1] - u[x,y+1] -2 *u[x+1,y] +2*u[x-1,y] -u[x+1, y+1] +u[x-1,y-1])*h/float(6)
            
            


            
            
    return newu



def G1stencil_T(u,h):
    
    """This routine builds the stencil for G1 matrix"""
    
    newu  = np.zeros((u.shape[0], u.shape[0]))
    
    for x in range(1, u.shape[0]-1):
        
        for y in range(1, u.shape[0]-1):
            

            newu[x,y] = (-u[x, y] +u[x-1,y])*h
            
            
    return newu


def G2stencil_T(u,h):
    
    """This rountine builds te stencil for G2 matrix"""
    
    newu  = np.zeros((u.shape[0], u.shape[0]))
    
    for x in range(1, u.shape[0]-1):
        
        for y in range(1, u.shape[0]-1):
            

        
            newu[x,y] = (-u[x,y]+u[x,y-1])*h
            


            
    return newu



            
    
def Uzawa(u,rhs,alpha, Ama, Lma):
    
    omega = 630
    
    
    (d,xdim,ydim) = u.shape
    
    h = 1/ float(xdim+2-1)
    
    newu  = np.zeros((d, xdim, ydim))
    
#    f1 = np.reshape(rhs[0], (xdim**2,))
#    
#    f2 = np.reshape(rhs[1], (xdim**2,))
#    
#    f3 = np.reshape(rhs[2], (xdim**2,))
    
    diff1 = rhs[0]- np.sqrt(alpha)* Lstencil(u[3])
    
    diff1 = np.reshape(diff1[1:-1,1:-1], ((xdim-2)**2,))
    
    diff2 = rhs[1] + G1stencil_T(u[3],h)
    
    diff2 = np.reshape(diff2[1:-1,1:-1], ((xdim-2)**2,))
    
    diff3 = rhs[2] + G2stencil_T(u[3],h)
    
    diff3 = np.reshape(diff3[1:-1,1:-1], ((xdim-2)**2,))
    
    Af = np.reshape(spsolve(Ama,diff1), ((xdim-2), (ydim-2)))
    
    Lf2 = np.reshape(spsolve(Lma, diff2), ((xdim-2),(ydim-2)))
    
    Lf3 = np.reshape(spsolve(Lma, diff3), ((xdim-2), (ydim-2)))
    
    newu[0][1:-1,1:-1] = Af
    
    newu[1][1:-1,1:-1] = Lf2
    
    newu[2][1:-1,1:-1] = Lf3 
    

        
    newu[3] = u[3] + omega*(np.sqrt(alpha)*Lstencil(newu[0]) - G1stencil(newu[1],h) - G2stencil(newu[2],h) -rhs[3])
    
    
    return newu


    