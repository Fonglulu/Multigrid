# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

The test currently builds a uniform grid of a given size, include the fem 
matrix. It then partitions the grid into a given number of sub-grids. The
sub-grids must contain all of the necessary information that will allow them
to be worked on in parallel. That includes, for example, the ghost node
table and neighbour node tables. Each sub grid is sent to a different
worker where the fem grid is plotted.

This file contains the worker code. It simply receives the subgrid
from the host processor and plots it.
"""

# Import appropriate information from other classes
from Communication import Communication
from PlotSquare import plot_fem_grid

    
#######################################################################
# cell
#
# This routine receives a subgrid from the host and plots it
#
#
# Input: Communication commun
#
# Output: Once the subgrids have been sent the routine returns
#
####################################################################### 
def worker(commun):
    grid = commun.get_grid(commun.get_host_no())
    #print '.......'
    plot_fem_grid(commun.get_my_no(), grid)
