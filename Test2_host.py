#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 20:02:35 2018

@author: shilu
"""

from grid.Grid import Grid
from BuildSquare import build_square_grid_matrix
from grid.function.FunctionStore import sin_soln, sin_rhs
from BuildEquation import build_equation_linear_2D
from BuildEquation import Poisson_tri_integrate
from grid.GridDivide import subdivide_grid
from Communication import Communication
from PlotSquare import plot_fem_grid

def host(commun):
    """Define the host processor"""

    # Initialise the grid
    grid = Grid()

    # Build a 5*5 grid
    i = 6
    n = 2**i+1
    
    # Specify the boundary conditions and the rhs
    true_soln = sin_soln
    rhs = sin_rhs
    

    
    # Find the number of processors
    p = commun.get_no_workers()

    # Specify the boundary conditions and the rhs
 
    # Create a square FEM grid of size n*n
    build_square_grid_matrix(n, grid, true_soln)

    # Calculate and store the A matrix and rhs vector
    build_equation_linear_2D(grid, Poisson_tri_integrate, rhs)

    # Subdivide the grids
    sub_grids = subdivide_grid(p, grid)
    
#    for cell in sub_grids:
#        plot_fem_grid(sub_grids[cell][0])

    # Send each subgrid to the specified worker processor
    for i in range(p):
        commun.send_grid(i, sub_grids[i][0])
 
    