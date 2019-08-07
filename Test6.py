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

"""

# Import appropriate information from other classes
from Communication import Communication
from Test6_worker import worker
from Test6_host import host

    
# put in class and set tags in class as well as comm
# and cell id for host

commun = Communication()   
if commun.get_my_no() == commun.get_host_no():
    host(commun)
else:
    worker(commun)