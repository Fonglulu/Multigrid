# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This class stores a table of FEM edges. See the Edge class for 
a description of a FEM edge. Access to the table is through the local 
node ids of the nodes sitting at the two endpoints of the edge.

The set and get operations make copies the information.

The iterator class, endpt_iterator, returns a reference to the edge.

"""

# Import appropriate information from other classes
from EdgeStar import EdgeStar
from copy import deepcopy


#######################################################################
# endpt_iterator
#
# Define an iterator generator that loops over all of the nodes joined
# by an edge to node_id.
#
# If there are no edges joined to node_id, then the routine will raise
# an error
#
# Input: EdgeTable edge_table
#        NodeID node_id
#
# Output: NodeID endpoint
#
#######################################################################
def endpt_iterator(edge_table, node_id):
    """Iterate over the end points of edges joined to node_id""" 
    
    # If there are no edges joined to the node then return
    if edge_table.get_no_endpoints(node_id) == 0:
        return
        
    # Find the edges star
    edge_star = edge_table._edge_stars[str(node_id)]

    # Loop over the end points
    for key, endpoint in edge_star._edge_ends.iteritems():
        yield endpoint.get_id()
        
        
###########################################################################   
class EdgeTable:
    """ A table of FEM edge"""
    
    # A dictionary (table) of edge stars
    _edge_stars = {}
    
    def __init__(self, edge_stars = dict()):
        """Initialise an edge table
        
        By default, creates an empty table
        """
        self._edge_stars = deepcopy(edge_stars)
     
    
    def set_boundary_function(self, endpt1, endpt2, boundary_func):
        """Set the boundary function"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Set the boundary function
        self._edge_stars[str(endpt1)].set_boundary_function(endpt2, 
                         boundary_func)
        
        
    def set_refine_type(self, endpt1, endpt2, refine_type):
        """Set the refinement type"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Set the refine type
        self._edge_stars[str(endpt1)].set_refine_type(endpt2, refine_type)
        
        
    def set_location(self, endpt1, endpt2, location):
        """Set where the edge sits (interior or boundary)"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Set the location
        self._edge_stars[str(endpt1)].set_location(endpt2, location)


    def set_error_indicator(self, endpt1, endpt2, indicator):
        """Set the error indicator"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Set the error indicator
        self._edge_stars[str(endpt1)].set_error_indicator(endpt2, indicator)

        
    def set_refine_level(self, endpt1, endpt2, level):
        """Set the refinement level when the edge was included in the domain"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Set the refine level
        self._edge_stars[str(endpt1)].set_refine_level(endpt2, level)


    def add_edge_all(self, edge):
        """Add a new edge to the table
        
        The edge is copied across
        
        """
        endpt1 = edge.get_endpt1()
        endpt2 = edge.get_endpt2()
        
        # Make sure the edge is not already in the table
        assert not self.is_edge(endpt1, endpt2), \
            "edge is already in table :"+str(endpt1)+"_"+str(endpt2)
            
        # If the edge star is not already in the table, create a new one
        if not self._edge_stars.has_key(str(endpt1)):
            self._edge_stars[str(endpt1)] = EdgeStar(endpt1)
            
        # Store the edge in an edge star
        self._edge_stars[str(endpt1)].add_endpoint(endpt2, 
                edge.get_location(), edge.get_boundary_function(), 
                edge.get_refine_type(), 
                edge.get_error_indicator(), edge.get_refine_level())


    def add_edge(self, endpt1, endpt2):
        """Add a new edge to the table
        
        The edge is copied across
        
        """
        
        # Make sure the edge is not already in the table
        assert not self.is_edge(endpt1, endpt2), \
            "edge is already in table :"+str(endpt1)+"_"+str(endpt2)
            
        # If the edge star is not already in the table, create a new one
        if not self._edge_stars.has_key(str(endpt1)):
            self._edge_stars[str(endpt1)] = EdgeStar(endpt1)
            
        # Store the edge in an edge star
        self._edge_stars[str(endpt1)].add_endpoint(endpt2)
            
            
    def delete_edge(self, endpt1, endpt2):
        """Delete an edge from the table"""

        # Make sure the edge is in the table
        assert self.is_edge(endpt1, endpt2), \
            "edge is node in table:"+str(endpt1)+"_"+str(endpt2)
        
        # Remove the edge from the dictionary
        self._edge_stars[str(endpt1)].delete_endpoint(endpt2)
        
        # If no more edges are joined to endpt1, remove the star
        if self._edge_stars[str(endpt1)].get_no_endpoints() == 0:
            del self._edge_stars[str(endpt1)]
        
       
    def is_edge(self, endpt1, endpt2):
        """Is the edge in the table?"""
        if len(self._edge_stars) == 0:
            return False
        if not self._edge_stars.has_key(str(endpt1)):
            return False
        return self._edge_stars[str(endpt1)].is_in(endpt2)


    def get_boundary_function(self, endpt1, endpt2):
        """Get the reference to the boundary function"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Get the boundary function
        return self._edge_stars[str(endpt1)].get_boundary_function(endpt2)
        
        
    def get_refine_type(self, endpt1, endpt2):
        """Get the refinement type"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Get the refinement type
        return self._edge_stars[str(endpt1)].get_refine_type(endpt2)
        
        
    def get_location(self, endpt1, endpt2):
        """Get the position where the edge sits in the domain"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Get the location
        return self._edge_stars[str(endpt1)].get_location(endpt2)


    def get_error_indicator(self, endpt1, endpt2):
        """Get the error indicator"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Get the error indicator
        return self._edge_stars[str(endpt1)].get_error_indicator(endpt2)


    def get_refine_level(self, endpt1, endpt2):
        """Get the refinement level when the edge was included in the domain"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Get the refinement level
        return self._edge_stars[str(endpt1)].get_refine_level(endpt2)
        
    
    def get_no_endpoints(self, endpt1):
        """Get the number of edges whose first node id id endpt1"""
        
        # If nothing is joined to the end point return 0
        if not self._edge_stars.has_key(str(endpt1)):
            return 0
            
        # Count the number of end points in the star
        return self._edge_stars[str(endpt1)].get_no_endpoints()


    def get_no_edges(self):
        """Get the total number of edges in the domain"""
        
        # Loop through and count the number of endpoints for each star
        no_edges = 0
        for node_id, star in self._edge_stars.iteritems():
            no_edges += star.get_no_endpoints()
        return no_edges
        
        
    def get_edge(self, endpt1, endpt2):
        """Extract the edge between ids endpt1 and endpt2 from the table"""
        
        # Make sure the edge is in the table
        assert self._edge_stars.has_key(str(endpt1)), \
            "edge not in table :"+str(endpt1)+"_"+str(endpt2)
            
        # Get the edge
        return self._edge_stars[str(endpt1)].get_edge(endpt2)
        
        
    def __iter__(self):
        """Set up an iterator class to loop over the edge stars in the table"""
        for staritem in self._edge_stars.iteritems():
            yield staritem[1]
            
    def display_edge_table(self):
        """Print the edge table
        
        The format is to loop through the table and print each edge
        
        """
        for star in self._edge_stars:
            self._edge_stars[star].display_star()
