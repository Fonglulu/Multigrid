#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 07:15:17 2018

@author: shilu
"""

# Import appropriate information from other classes
from grid.Grid import Grid
from BuildSquare import build_square_grid
from PlotSquare import plot_fem_grid
from grid.function.FunctionStore import exp_soln,exp_rhs, sin_soln, sin_rhs
from BuildEquation import build_equation_linear_2D,\
    Poisson_tri_integrate
from grid.GridDivide import subdivide_grid


# Initialise the grid
grid = Grid()




# Build a 5*5 grid
i = 3
n = 2**i + 1

# Specify the boundary condition and the rhs
true_soln = sin_soln
rhs = sin_rhs

# Create a square FEM grid of size n*n
build_square_grid(n, grid, true_soln)

# Calculate and store the A matrix and rhs vector
build_equation_linear_2D(grid, Poisson_tri_integrate, rhs)

# Need to add neighbour nodes and ghost grids
sub_grids = subdivide_grid(2,grid)

# Plot the solution
for cell in sub_grids:
    plot_fem_grid(cell, sub_grids[cell][0])