#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 16:57:14 2018

@author: shilu
"""

from numpy import np

from Rich_uniform import Astencil, G1stencil, G2stencil, Lstencil


def GS(u, rhs,alpha):
    
    pass
    
    
    [xdim,ydim] = rhs[0][1:-1, 1:-1].shape
    
    h = 1/ float(xdim+2-1)
    
    
    
    crhs = rhs[0] + G1stencil(h, u[1]) + G2stencil(h, u[2])
    
    
    
    