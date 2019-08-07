# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

The test currently builds a uniform grid of a given size, include the fem 
matrix. It then partitions the grid into a given number of sub-grids. The
sub-grids must contain all of the necessary information that will allow them
to be worked on in parallel. That includes, for example, the ghost node
table and neighbour node tables. Each sub grid is sent to a different
worker)where the fem grid is plotted.

This file contains the host or master code. This is where the original 
grid is built and the subgrids are created. Each subgrid is then sent to
a worker.

"""

# Import appropriate information from other modules
from grid.Grid import Grid
from BuildSquare import build_square_grid_matrix
from grid.function.FunctionStore import exp_soln, exp_rhs, sin_soln, sin_rhs
from BuildEquation import build_equation_linear_2D
from BuildEquation import Poisson_tri_integrate
from grid.GridDivide import subdivide_grid
from Communication import Communication
    

#######################################################################
# host
#
# This routine firstly builds the original full domain grid. It then
# used metis to subdivie the grid. Each subgrid is sent to a worker
# processor
#
#
# Input: Communication commun
#
# Output: Once the subgrids have been sent the routine returns
#
####################################################################### 
def host(commun):
    """Define the host processor"""

    # Initialise the grid
    grid = Grid()

    # Build a 5*5 grid
    i = 4
    n = 2**i+1
    
    # Find the number of processors
    p = commun.get_no_workers()

    # Specify the boundary conditions and the rhs
    true_soln = sin_soln
    rhs = sin_rhs

    # Create a square FEM grid of size n*n
    build_square_grid_matrix(n, grid, true_soln)

    # Calculate and store the A matrix and rhs vector
    build_equation_linear_2D(grid, Poisson_tri_integrate, rhs)

    # Subdivide the grids
    sub_grids = subdivide_grid(p, grid)

    # Send each subgrid to the specified worker processor
    for i in range(p):
        commun.send_grid(i, sub_grids[i][0])
        print '...........'

