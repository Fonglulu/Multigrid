# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This class stores a table of FEM algebraic connection. 

"""

# Import appropriate information from other classes
from ConnectStar import ConnectStar
from copy import deepcopy


#######################################################################
# connect_iterator
#
# Define an iterator generator that loops over all of the nodes joined
# by an connection to node_id.
#
# If there are no connections to node_id, then the routine will raise
# an error
#
# Input: ConnectTable connection_table
#        NodeID node_id
#
# Output: NodeID endpoint
#
#######################################################################
def connect_iterator(connection_table, node_id):
    """Iterate over the end points of connections joined to node_id""" 
    
    # If there are no connections joined to the node then return
    if connection_table.get_no_endpoints(node_id) == 0:
        return

    # Find the star joined to node_id
    connection_star = connection_table._connection_stars[str(node_id)]

    # Loop over the connections joined to node_id
    for key, endpoint in connection_star._connection_ends.iteritems():
        yield endpoint.get_id()
        
        
###########################################################################   
class ConnectTable:
    """ A table of FEM algebraic connection"""
    
    # A dictionary (table) of connection stars
    _connection_stars = {}
    
    def __init__(self, connection_stars = dict()):
        """Initialise a connection table
        
        By default, creates an empty table
        """
        self._connection_stars = deepcopy(connection_stars)
     
    
    def set_matrix_value(self, endpt1, endpt2, value):
        """Set the value of the corresponding entry in the matrix"""

        # Make sure the connection is in the table
        assert self._connection_stars.has_key(str(endpt1)), \
            "connection not in table :"+str(endpt1)+"_"+str(endpt2)

        # Set the value
        self._connection_stars[str(endpt1)].set_value(endpt2, value)
        

    def add_connection(self, endpt1, endpt2, value= 0.0):
        """Add a new connection to the table
        
        The connection is copied across
        
        """
        
        # Make sure the connection is not already in the table
        assert not self.is_connection(endpt1, endpt2), \
            "connection is already in table :"+str(endpt1)+"_"+str(endpt2)
        
        # If the star is not in the table, create a new star
        if not self._connection_stars.has_key(str(endpt1)):
            self._connection_stars[str(endpt1)] = ConnectStar(endpt1)
            
        # Store the connection in an connection star
        self._connection_stars[str(endpt1)].add_connection(endpt2, value)
        
        
    def delete_connection(self, endpt1, endpt2):
        """Delete an connection from the table"""

        # Make sure the connection is in the table
        assert self.is_connection(endpt1, endpt2), \
            "connection is node in table:"+str(endpt1)+"_"+str(endpt2)
        
        # Remove the connection from the dictionary
        self._connection_stars[str(endpt1)].delete_connection(endpt2)
        
        # If no more connections are joined to endpt1, remove the star
        if self._connection_stars[str(endpt1)].get_no_endpoints() == 0:
            del self._connection_stars[str(endpt1)]
        
       
    def is_connection(self, endpt1, endpt2):
        """Is the connection in the table?"""

        # If the table is empty, then return False
        if len(self._connection_stars) == 0:
            return False

        # If there are no stars joined to endpt1, then return False
        if not self._connection_stars.has_key(str(endpt1)):
            return False

        # Check if endpt2 is connected to endpt1
        return self._connection_stars[str(endpt1)].is_in(endpt2)


    def get_matrix_value(self, endpt1, endpt2):
        """Get the corresponding value in the matrix"""

        # Make sure the connection is in the table
        assert self._connection_stars.has_key(str(endpt1)), \
            "connection not in table :"+str(endpt1)+"_"+str(endpt2)

        # Return the value
        return self._connection_stars[str(endpt1)].get_value(endpt2)
        
        
    def get_no_endpoints(self, endpt1):
        """Get the number of connections whos first node id id endpt1"""

        # Make sure the connection is in the table
        assert self._connection_stars.has_key(str(endpt1)), \
            "connection star not in table :"+str(endpt1)

        # Return the number of endpoints
        return self._connection_stars[str(endpt1)].get_no_endpoints()


    def get_no_connections(self):
        """Get the total number of connections in the domain"""

        # Initialise the number of connection to 0
        no_connections = 0

        # Loop through the stars and add in the number of connections to each star
        for node_id, star in self._connection_stars.iteritems():
            no_connections += star.get_no_endpoints()
        return no_connections
          
        
    def __iter__(self):
        """Set up an iterator class to loop over the stars in the table"""
        for staritem in self._connection_stars.iteritems():
            yield staritem[1]
            
    def display_connection_table(self):
        """Print the algebraic connection table table
        
        The format is to loop through the table and print each connection
        
        """
        for star in self._connection_stars:
            self._connection_stars[star].display_star()
