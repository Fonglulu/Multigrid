# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 22:20:27 2013

@author: stals

These classes keep a record of a FEM edge stars. An edge star is a list of
all of the edges joined to a given node. For a description of the information
stored with each edge in the star see the Edge class.

The set and get operations copies the information across, eg 
get_endpt_id() returns a local copy of the id and 
set_endpt_id() copies the information across (uses python copy command).
"""

# Import appropriate information from other classes
from NodeID import NodeID
from Edge import DomainSet, RefineSet, Edge
from function.FunctionStore import FunctionStore, zero
from copy import copy, deepcopy
   

###########################################################################   
class EdgeEnd:
    """ A finite element edge"""
    
    # The id of the second node in the edge
    _node_id = NodeID()
    
    # Where the edge sits in the domain
    _location = DomainSet.interior
    
    # The boundary function (dummy values for the interior edges)
    _boundary_function = zero
    
    # What is the edge type in terms of the refinement routine
    _refine_type = RefineSet.not_base_edge
    
    # What is the error indicator for the given edge
    _error_indicator = 0.0
    
    # This edge was created at what level of refinement
    _refine_level = 0
    
   
    def __init__(self, node_id = NodeID(), location = DomainSet.interior, 
            boundary_func = zero, refine_type = RefineSet.not_base_edge, 
            error_indicator = 0.0, refine_level = 0):
        """ Initialise the edge
        
        By default, a dummy node id is assigned to the endpoint of the
        edge. It is assumed to be an interior edge and is not a base type edge.
        The error indicator is 0.0 and the refinement level is 0
        
        """
        self._node_id = copy(node_id)
        self._location = location
        self._boundary_function = boundary_func
        self._refine_type = refine_type
        self._error_indicator = error_indicator
        self._refine_level = refine_level       
        
    def set_id(self, node_id):
        """Set the id of the second end node of the edge"""
        self._node_id = copy(node_id)
        
    def set_location(self, location):
        """Set where the function sits (interior or boundary)"""
        self._location = location
        
    def set_error_indicator(self, indicator):
        """Set the error indicator"""
        self._error_indicator = indicator
        
    def set_refine_level(self, level):
        """Set the refinement level when the edge was included in the domain"""
        self._refine_level = level
        
    def set_boundary_function(self, boundary_func):
        """Set the boundary function"""
        self._boundary_function = boundary_func
        
    def set_refine_type(self, refine_type):
        """Set the refinement type"""
        self._refine_type = refine_type
        
    def get_id(self):
        """Get the id second end node of the edge"""
        node_id = copy(self._node_id)
        return node_id
        
    def get_location(self):
        """Get the position where the edge sits in the domain"""
        return self._location
        
    def get_error_indicator(self):
        """Get the error indicator"""
        return self._error_indicator
        
    def get_refine_level(self):
        """Get the refinement level when the edge was included in the domain"""
        return self._refine_level
        
    def get_boundary_function(self):
        """Get the reference to the boundary function"""
        return self._boundary_function
        
    def get_refine_type(self):
        """Get the refinement type"""
        return self._refine_type
       

########################################################################### 
class EdgeStar:
    """A collection of edges joined to a given node"""
    
    # The id of the node at the centre of the star
    _node_id = NodeID()
    
    # A collection of edges joined to the node
    _edge_ends = {}
    

    def __init__(self, node_id=NodeID(), edge_ends = dict()):
        """ Initialise the star
         
        By default, the centre node is assigned a dummy id and the set of
        edges is empty
        
        """
        self._node_id = copy(node_id)
        self._edge_ends = deepcopy(edge_ends) 
        
    def set_id(self, node_id):
        """Set the id of the node at the centre of the star"""
        self._node_id = copy(node_id)
        
    def set_boundary_function(self, endpt, boundary_func):
        """Set the boundary function"""
        # A copy of the function store
        func_store = FunctionStore()
        assert func_store.is_in_store(boundary_func), \
            "the boundary function must be in FunctionStore "
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        self._edge_ends[str(endpt)].set_boundary_function(boundary_func)
        
    def set_refine_type(self, endpt, refine_type):
        """Set refinement type"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        self._edge_ends[str(endpt)].set_refine_type(refine_type)
        
    def set_location(self, endpt, location):
        """Set where the edge sits (interior or boundary)"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        self._edge_ends[str(endpt)].set_location(location)

        
    def set_error_indicator(self, endpt, indicator):
        """Set the error indicator"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        self._edge_ends[str(endpt)].set_error_indicator(indicator)

        
    def set_refine_level(self, endpt, level):
        """Set the refinement level when the edge was included in the domain"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        self._edge_ends[str(endpt)].set_refine_level(level)

    def add_endpoint(self, endpt, location = DomainSet.interior, \
            boundary_func = zero, refine_type = RefineSet.not_base_edge, \
            error_indicator = 0.0, refine_level = 0):
        """Add an additional edge to the collection of edges"""
        edge_end = EdgeEnd(endpt, location, boundary_func, 
                           refine_type, error_indicator, refine_level)
        self._edge_ends[str(endpt)] = copy(edge_end)
        
    def delete_endpoint(self, endpt):
        """Delete an edge from the star"""
        assert self.is_in(endpt), \
            "end point not in table :"+str(endpt)
        del self._edge_ends[str(endpt)]
        
    def is_in(self, endpt):
        """Is the edge in the collection of edges?"""
        if (len(self._edge_ends) == 0):
            return False
        return self._edge_ends.has_key(str(endpt))
        
    def get_id(self):
        """Get the id of the node at the centre of the star"""
        node_id = copy(self._node_id)
        return node_id
        
    def get_boundary_function(self, endpt):
        """Get the boundary function"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        return self._edge_ends[str(endpt)].get_boundary_function()
        
    def get_refine_type(self, endpt):
        """Get the refinement type"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        return self._edge_ends[str(endpt)].get_refine_type()
        
    def get_location(self, endpt):
        """Get the position where the edge sits in the domain"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        return self._edge_ends[str(endpt)].get_location()
        
    def get_error_indicator(self, endpt):
        """Get the error indicator"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        return self._edge_ends[str(endpt)].get_error_indicator()
        
    def get_refine_level(self, endpt):
        """Get the refinement level when the edge was included in the domain"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        return self._edge_ends[str(endpt)].get_refine_level()

    def get_no_endpoints(self):
        """How many edges are connected to the centre of the star?"""
        return len(self._edge_ends)
        
    def get_edge(self, endpt):
        """Get the edge going from the centre of the star to the given endpt"""
        assert self._edge_ends.has_key(str(endpt)), \
            "edge not in table :"+str(self._node_id)+"_"+str(endpt)
        end_point = self._edge_ends[str(endpt)]
        edge = Edge(self._node_id, endpt, end_point.get_location(), 
                end_point.get_boundary_function(), 
                end_point.get_refine_type(), 
                end_point.get_error_indicator(), 
                end_point.get_refine_level())
        return edge
        
    def __iter__(self):
        """Set up an iterator class to loop over the edges in the star"""
        for enditem in self._edge_ends.iteritems():
            yield enditem[1]
            
    def display_star(self):
        """Print the edge star
        
        The format is to loop through the edges in the star and print both 
        endpoints
        
        """
        for end in self:
            print str(self._node_id) + "_" + str(end.get_id())
        