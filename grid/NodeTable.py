# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This class stores a table of FEM nodes. See the Node class for 
a description of a FEM node. Access to the table is generally
through the local node id, although there it is possible to 
map between the local and global ids.

The set and get operations copies the information across, eg 
get_node_id() returns a local copy of the id and 
set_node_id() copies the information across (uses python copy 
command).

The iterator class, node_iterator, returns a reference to the node.
"""

# Import appropriate information from other classes

from NodeID import NodeID
from Node import Node
from copy import copy, deepcopy


#######################################################################
# NodeIterator
#
# Define an iterator generator that loops over the nodes in a 
# NodeTable
#
# Input: NodeTable node_table
#
# Output: Node node
#
#######################################################################
def node_iterator(node_table):
    """Loop over the nodes in the node table"""
    for key, node in node_table._nodes.iteritems():
        yield node
 


def not_slave_node(grid):
    """ Loop over the interior nodes in the node table"""
    
    for node in node_iterator(grid):
        if not node.get_slave():
            yield node
###########################################################################   
       
class NodeTable:
    """ A table of FEM nodes"""
    
    # Node dimension (should match up with the dimension of the
    # nodes)
    _dim = 0
    
    # A dictionary (table) of nodes
    _nodes = {}
    
    # A dictionary used to map a global id into a local id
    _id_map = {}
    
    # A static variable keep track of the number of nodes that
    # have been created
    _counter = 0 
    
    def __init__(self, dim=2, nodes = dict(), id_map = dict()):
        """Initialise a node table
        
        By default, creates an empty two dimensional table
        """
        self._dim = dim
        self._nodes = deepcopy(nodes)
        self._id_map = deepcopy(id_map)
        
    def set_slave(self, node_id, is_slave):
        """Is the node a slave node (T/F)"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Set the condition of the node
        self._nodes[str(node_id)].set_slave(is_slave)

    def set_value(self, node_id, value):
        """Set current value"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id) 
            
        # Set the node value
        self._nodes[str(node_id)].set_value(value)
        
    def set_coord(self, node_id, coord):
        """Set node coordinates"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Set the node coordinates
        self._nodes[str(node_id)].set_coord(coord)

    def set_load(self, node_id, load):
        """Set current load"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Set the load (or RHS) value
        self._nodes[str(node_id)].set_load(load)

        
    def create_node(self, global_id, coord, is_slave, value):
        """Create a new node
        
        This routine should be used to create new nodes that will
        be added to the node table. The routine is responsible for
        assigning the local id, so if nodes are created outside of
        this routine care needs to be take with regards to the local
        id. It is assumed that the caller has determined the global
        id (probably based on the node coordinates). The dimension
        of the node is the same as the node table
        
        """
        # Initialise the node
        node = Node(self._dim)
        node_id = NodeID()
        
        # The local id is just a counter
        node_id.set_no(NodeTable._counter)
        NodeTable._counter = NodeTable._counter+1
        node.set_node_id(node_id)
        
        # Copy across the global id
        node.set_global_id(global_id)
        
        # Set the coordinates, checking that dimension is correct
        assert self._dim == len(coord), \
            "coordinates should be of length %d"% self._dim 
        node.set_coord(coord)
        
        # Record if the node is a slave node
        node.set_slave(is_slave)
        
        # Record the node value
        node.set_value(value)
        
        return node
        
    def add_node(self, node):
        """Add a new node to the table
        
        The node is copied across
        
        """
        # Make sure the node is not already in the table
        assert not self.is_in(node.get_node_id()), \
            "node is already in table :"+str(node.get_node_id())
            
        # Store the node in a dictionary using the local id as
        # a key
        self._nodes[str(node.get_node_id())] = copy(node)
        
        # Use a second dictionary to map from the global id to
        # the local id
        self._id_map[str(node.get_global_id())] \
            = copy(node.get_node_id())
        
    def delete_node(self, node_id):
        """Delete a node from the table"""
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Give the local id find the global id
        global_id = self.get_global_id(node_id)
        
        # Remove the node from the dictionary
        del self._nodes[str(node_id)]
        
        # Remove the node from the id map
        del self._id_map[str(global_id)]
       
    def is_in(self, node_id):
        """Is the node in the table?"""
        
        # If the table is empty return False
        if len(self._nodes) == 0:
            return False
            
        # Else check to see if the node is in the dictionary
        return self._nodes.has_key(str(node_id))
        
    def is_in_global(self, global_node_id):
        """Is the node in the table?"""
        
        # If the table is empty return False
        if len(self._id_map) == 0:
            return False
            
        # Else check to see if the node is in the dictionary
        return self._id_map.has_key(str(global_node_id))

    def get_no_nodes(self):
        """Get the number of nodes in the table"""
        return len(self._nodes)
        
    def get_dim(self):
        """Get the dimension of the FEM domain"""
        return self._dim
        
    def get_slave(self, node_id):
        """Is the node a slave node?"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Get the node condition
        return self._nodes[str(node_id)].get_slave()

    def get_value(self, node_id):
        """Get the current node value"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Get the node value
        return self._nodes[str(node_id)].get_value()

    def get_load(self, node_id):
        """Get the current load value"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Get the load (or RHS) value
        return self._nodes[str(node_id)].get_load()  
     
    def get_coord(self, node_id):
        """Get the node coordinates"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Get the node's coordinates
        return self._nodes[str(node_id)].get_coord()   

    def get_global_id(self, node_id):
        """Given a local id, what is the corresponding global id"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Get the global node id
        return self._nodes[str(node_id)].get_global_id()
        
    def get_local_id(self, global_id):
        """Given a global id, what is the corresponding local id"""
        
        # Make sure the node is in the table
        assert self._id_map.has_key(str(global_id)), \
            "node not in table :"+str(global_id)   
            
        # Get the local id
        return self._id_map[str(global_id)]
        
    def get_node(self, node_id):
        """Get the FEM node with local id"""
        
        # Make sure the node is in the table
        assert self.is_in(node_id), \
            "node not in table :"+str(node_id)
            
        # Return a copy of the node
        node = copy(self._nodes[str(node_id)])
        return node

    def __iter__(self):
        """Set up an iterator class to loop over the nodes in the table"""
        for nodeitem in self._nodes.iteritems():
            yield nodeitem[1]
            
    def display_node_table(self):
        """Print the node table
        
        The format is to loop through the table and print each node
        
        """
        for node in self:
            print(node)

    