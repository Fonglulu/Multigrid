# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:35:04 2013

@author: stals

This class keeps a record a FEM node and provides some elementary 
operations on the node

The set and get operations copies the information across, eg 
get_node_id() returns a local copy of the id and 
set_node_id() copies the information across (uses python copy 
command).
"""

# Import appropriate information from other classes

from NodeID import NodeID
from GlobalNodeID import GlobalNodeID
from numpy import array
from copy import copy

class Node:
    """ A finte element node"""
    
    # Node dimension (should match up with the dimension of the
    # node table)
    _dim = 0
    
    # Local node id
    _node_id = NodeID()
    
    # Global node id
    _global_id = GlobalNodeID()
    
    # Is the node a slave node (i.e. does it sit on the boundary?)
    _slave = False
    
    # What is the value assigned to the node
    _value = 0.0
    
    # What is the load (RHS value)
    _load = 0.0
    
    # What are the node's coordinates (size should be same as dim)
    _coord = array([])
   
   
    def __init__(self, dim=2, coord=[0, 0], 
                 slave=False, value=0.0, load=0.0):
        """Initialise a node
        
        By default, the node is a two dimensional node with 
        coordinates (0,0), value 0 and load 0. It is also assumed
        to be an interior node (i.e. not a slave node)
        """
        self._dim = dim
        self._coord = coord
        self._slave = slave
        self._value = value
        self._load = load
        self._node_id = NodeID()
        self._global_id = GlobalNodeID()
        
    def __str__(self):
        """Convert a global node id into a string so it can be printed
        
        The format is the local id (global id) value
        
        """
        return str(self._node_id) + \
            "  (" + str(self._global_id) + ") "+ \
            repr(self._value)
    
     
    def set_node_id(self, node_id):
        """Set the local node id
        
        The local id is copied across
        
        """
        self._node_id = copy(node_id)
    
    def set_global_id(self, global_id):
        """Set the global node id
        
        The global id is copied across
        
        """
        self._global_id = copy(global_id)
        
    def set_dim(self, dim):
        """Set the node dimension"""
        self._dim = dim
        
    def set_load(self, load):
        """Set the load (RHS value)"""
        self._load = load
        
    def set_value(self, value):
        """Set current value"""
        self._value = value
    
    def set_slave(self, slave):
        """Is the node a slave node (T/F)"""
        self._slave = slave
        
    def set_coord(self, coord):
        """Set the coordinates
        
        The length of the coordinates should agree with the
        node dimension
        
        """
        assert self._dim == len(coord), \
            "coordinates should be of length %d"% self._dim 
        self._coord = coord[:]

        
    def get_node_id(self):
        """Get the local node id"""
        node_id = copy(self._node_id)
        return node_id
        
    def get_global_id(self):
        """Get the global node id"""
        global_id = copy(self._global_id)
        return global_id 
    
    def get_dim(self):
        """Get the node dimension"""
        return self._dim 
       
    def get_load(self):
        """Get the load (RHS value)"""
        return self._load
        
    def get_value(self):
        """Get the node value"""
        return self._value
        
    def get_slave(self):
        """Is the node a slave node (T/F)"""
        return self._slave
    
    def get_coord(self):
        """Get the node coordinates"""
        coord = self._coord[:]
        return coord
        

        