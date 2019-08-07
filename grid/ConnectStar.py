# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 22:20:27 2013

@author: stals

These classes keep a record of the FEM connection stars. An connection star 
is a list of all of the algebraic connections to a given node. 

The set and get operations copies the information across, eg 
get_endpt_id() returns a local copy of the id and 
set_endpt_id() copies the information across (uses python copy command).
"""

# Import appropriate information from other classes
from NodeID import NodeID
from copy import copy, deepcopy
   

###########################################################################   
class ConnectEnd:
    """ A finite element connection"""
    
    # The id of the second node in the connection
    _node_id = NodeID()
    
    # The corresponding matrix value
    _value = 0.0
   
    def __init__(self, node_id = NodeID(), value = 0.0):
        """ Initialise the connection
        
        By default, a dummy node id is assigned to the endpoint of the
        connection and it is given a value of 0
        
        """
        self._node_id = copy(node_id)
        self._value = value 
        
    def set_id(self, node_id):
        """Set the id of the second end node of the connection"""
        self._node_id = copy(node_id)
        
    def set_value(self, value):
        """Set the corresponding matrix value"""
        self._value = value
        
    def get_id(self):
        """Get the id second end node of the connection"""
        node_id = copy(self._node_id)
        return node_id
        
    def get_value(self):
        """Get the corresponding matrix value"""
        return self._value
       

###########################################################################
class ConnectStar:
    """A collection of connections joined to a given node"""
    
    # The id of the node at the centre of the star
    _node_id = NodeID()
    
    # A collection of connections joined to the node
    _connection_ends = {}
    

    def __init__(self, node_id=NodeID(), connection_ends = dict()):
        """ Initialise the star
        
        By default, the centre node is assigned a dummy id and the set of
        connections is empty
        
        """
        self._node_id = copy(node_id)
        self._connection_ends = deepcopy(connection_ends) 
        
    def set_id(self, node_id):
        """Set the id of the node at the centre of the star"""
        self._node_id = copy(node_id)
        
    def set_value(self, endpt, value):
        """Set the corresponding matrix value"""
        assert self._connection_ends.has_key(str(endpt)), \
            "connection not in table :"+str(self._node_id)+"_"+str(endpt)
        self._connection_ends[str(endpt)].set_value(value)
        

    def add_connection(self, endpt, value = 0.0):
        """Add an additional connection to the collection of connections"""
        connection_end = ConnectEnd(endpt, value)
        self._connection_ends[str(endpt)] = copy(connection_end)
        
    def delete_connection(self, endpt):
        """Delete an connection from the star"""
        assert self.is_in(endpt), \
            "end point not in table :"+str(endpt)
        del self._connection_ends[str(endpt)]
        
    
    def is_in(self, endpt):
        """Is the connection in the collection of connections?"""
        if (len(self._connection_ends) == 0):
            return False
        return self._connection_ends.has_key(str(endpt))
        
    def get_id(self):
        """Get the id of the node at the centre of the star"""
        node_id = copy(self._node_id)
        return node_id
        
    def get_value(self, endpt):
        """Get the corresponding matrix value"""
        assert self._connection_ends.has_key(str(endpt)), \
            "connection not in table :"+str(self._node_id)+"_"+str(endpt)
        return self._connection_ends[str(endpt)].get_value()
        
    def get_no_endpoints(self):
        """How many connections are connected to the centre of the star?"""
        return len(self._connection_ends)
        
    def __iter__(self):
        """Set up an iterator class to loop over the connections in the star"""
        for enditem in self._connection_ends.iteritems():
            yield enditem[1]
            
    def display_star(self):
        """Print the connection star
        
        The format is to loop through the edges in the star and print both 
        endpoints
        
        """
        for end in self:
            print str(self._node_id) + "_" + str(end.get_id())
            
