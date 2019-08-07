# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This class stores a finite element grid. It inherits most of the
operations from the NodeTable, EdgeTable and ConnectTable classes.
The exceptions are the ghost node table and neighbour node table that
are stored as variables with the Grid class.

"""

# import appropriate information from other classes
from EdgeTable import EdgeTable
from NodeTable import NodeTable
from ConnectTable import ConnectTable
from NghNodes import NghNodes

class Grid(NodeTable, EdgeTable, ConnectTable):
    """ A finite element grid"""
    
    _ghost_table = NodeTable()
    _ghost_commun = NghNodes()
    _full_commun = NghNodes()
    
    def __init__(self, dim = 2):
        """Initialise a edge table
        
        By default, creates an empty table
        """
        NodeTable.__init__(self, dim)
        EdgeTable.__init__(self)
        ConnectTable.__init__(self)    
        self._ghost_table = NodeTable(dim)
        self._ghost_commun = NghNodes()
        self._full_commun = NghNodes()
        
    def reference_ghost_table(self):
        """Return a reference to the ghost node table"""
        return self._ghost_table
        
    def reference_ghost_commun(self):
        """Return a reference the cells containing the full node copy of 
        the ghost nodes"""
        return self._ghost_commun
        
    def reference_full_commun(self):
        """Return a reference the cells containing any ghost node copy of 
        the full nodes"""
        return self._full_commun
        
        