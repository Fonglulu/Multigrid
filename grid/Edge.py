# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 22:20:27 2013

@author: stals

This class keeps a record of a FEM edge and provides some
elementary operations. This class is more of an interface class,
a nicer way of getting information in and out of the EdgeTable class. The
EdgeTable class itself is built on the EdgeStar class.

Note that edges have directions. That is, the edge between node 1 and node 2 is
not, in general, the same as the edge between node 2 and node 1.

The set and get operations copies the information across, eg 
get_endpt_id() returns a local copy of the id and 
set_endpt_id() copies the information across (uses python copy command).

The code currently only implements Dirichlet boundary conditions, to handle
Neumann boundaries the DomainSet class should be extended to include 
Dirichlet and Neumann boundaries.

"""

# Import appropriate information from other classes
from NodeID import NodeID
from function.FunctionStore import FunctionStore, zero
from copy import copy

class DomainSet:
    """Record where the edge sits in the domain.
    
    Does the edge sit in the interior of the domain or on the boundary of
    the domain
    
    """
    interior = 1
    boundary = 2
   
class RefineSet:
    """What is the edge type in terms of the refinement algorithm.
    
    Is the edge a base edge (sits opposite the newest node), an interface
    base edge (sits between two different levels of refinement) or not a base
    edge (neither base_edge or interface_base_edge)
    
    """
    base_edge = 1
    not_base_edge = 2
    interface_base_edge = 3
   
class Edge:
    """ A finite element edge"""
    
    # The node id of one end of the edge
    _endpt1 = NodeID()
    
    # The node id of the other end of the edge
    _endpt2 = NodeID()
    
    # Where the edge sits in the domain (boundary or interior)
    _location = DomainSet.interior
    
    # The boundary function (dummy value for interior nodes)
    _boundary_function = zero
    
    # What is the edge type in terms of the refinement routine
    _refine_type = RefineSet.not_base_edge
    
    # What is the error indicator for the given edge
    _error_indicator = 0.0
    
    # This edge was created at what level of refinement
    _refine_level = 0
   
    def __init__(self, endpt1=NodeID(), endpt2=NodeID(), 
                 location = DomainSet.interior, boundary_func = zero, 
                 refine_type = RefineSet.not_base_edge, error_indicator= 0.0,
                 refine_level = 0):
        """ Initialise the edge
        
        By default, two dummy node ids are assigned to the endpoints of the
        edge. It is assumed to be an interior edge and is not a base type edge.
        The error indicator is 0.0 and the refinement level is 0
        
        """
        self._endpt1 = copy(endpt1)
        self._endpt2 = copy(endpt2)
        self._location = location
        self._boundary_function = boundary_func
        self._refine_type = refine_type
        self._error_indicator = error_indicator
        self._refine_level = refine_level

    def __str__(self):
        """Convert an edge into a string so it can be printed
        
        The format is the id of the first endpoint _ id of the second endpoint 
        
        """
        return  str(self._endpt1) + " " + str(self._endpt2)
        
    def set_endpt1(self, node_id):
        """Set the id of the first end node of the edge"""
        self._endpt1  = copy(node_id)
        
    def set_endpt2(self, node_id):
        """Set the id of the second end node of the edge"""
        self._endpt2 = copy(node_id)
        
    def set_boundary_function(self, boundary_func):
        """Set the boundary function"""
        # A copy of the function store
        func_store = FunctionStore()
        assert func_store.is_in_store(boundary_func), \
            "the boundary function must be in FunctionStore "
        self._boundary_function = boundary_func
        
    def set_refine_type(self, refine_type):
        """Set the refinement type"""
        self._refine_type = refine_type
        
    def set_location(self, location):
        """Set where the edge sits (interior or boundary)"""
        self._location = location
        
    def set_error_indicator(self, indicator):
        """Set the error indicator"""
        self._error_indicator = indicator
        
    def set_refine_level(self, level):
        """Set the refinement level when the edge was included in the domain"""
        self._refine_level = level
        
    def get_endpt1(self):
        """Get the id of the first end node of the edge"""
        node_id = copy(self._endpt1)
        return node_id
        
    def get_endpt2(self):
        """Get the id of the second end node of the edge"""
        node_id = copy(self._endpt2)
        return node_id
        
    def get_location(self):
        """Get the position where the edge sits in the domain"""
        return self._location
        
    def get_boundary_function(self):
        """Get the reference to the boundary function"""
        return self._boundary_function
        
    def get_refine_type(self):
        """Get the refinement type"""
        return self._refine_type
        
    def get_error_indicator(self):
        """Get the error indicator"""
        return self._error_indicator
        
    def get_refine_level(self):
        """Get the refinement level when the edge was included in the domain"""
        return self._refine_level
        
        
 
        
