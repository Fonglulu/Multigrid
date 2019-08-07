#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 19:54:57 2018

@author: shilu
"""
from PCG import conjugate_gradient
from Preconditioner import identity_preconditioner
from DiagonalPrecond import DD_preconditioner

def worker(commun):
    grid = commun.get_grid(commun.get_host_no())
    
    
    loops, resid = conjugate_gradient(commun, grid, DD_preconditioner, 
                                      tolerance = 1.0E-6)
    
    
    print "CG method took " + str(loops)+" to reach a tolerance of "\
    +str(resid)

