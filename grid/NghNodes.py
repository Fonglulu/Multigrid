# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:19:15 2013

@author: stals

This module defines the classes used to keep track of the communication 
pattern.

Each ghost node needs to know the position of the corresponding full node
and each full node needs to know the position of any corresponding ghost node.
That information is stored in the record defined here.
"""

# import appropriate information from other classes


from copy import deepcopy, copy


#######################################################################
# ngh_iterator
#
# Define an iterator generator that loops over the neighbour cells
#
# Input: NghNodes ngh_nodes
#
# Output: neighbour_no, list of node ids in neighbouring processor
#
#######################################################################
def ngh_iterator(ngh_nodes):
    """Loop over the nodes in the node table"""
    for key, id_list in ngh_nodes._ngh_ids.iteritems():
        yield key, id_list
 
###########################################################################   
       
class NghNodes:
    """ A table of nodes sitting in neighbouring processor"""
    
    # A dictionary (table) of node ids in neighbour processor
    _ngh_ids = {}
    
    
    def __init__(self, ngh_ids = dict()):
        """Initialise a ngh node table
        
        By default, creates an empty table
        """
        self._ngh_ids = deepcopy(ngh_ids)
        
    def add_neighbour_node(self, ngh_no, node_id):
        """Add the node id to the list of neighbour nodes in cell ngh_no"""
        
        if not self.is_neighbour(ngh_no):
            self._ngh_ids[str(ngh_no)] = list()
        self._ngh_ids[str(ngh_no)].append(copy(node_id))
        
    def delete_neighbour_node(self, ngh_no, node_id):
        """Remove the node id from the list of neighbour nodes in cell ngh_no"""
        
        assert self.is_neighbour_node(ngh_no, node_id), \
            "node not in table :"+str(node_id)
        ngh_str = str(ngh_no)
        index = self._ngh_ids[ngh_str].index(node_id)
        self._ngh_ids[ngh_str].pop(index)
        if len(self._ngh_ids[ngh_str]) == 0:
            del self._ngh_ids[ngh_str]

    def is_neighbour(self, ngh_no):
        """Is the processor a neighbouring process?"""
        if len(self._ngh_ids) == 0:
            return False
        return self._ngh_ids.has_key(str(ngh_no))
        
    def is_neighbour_node(self, ngh_no, node_id):
        """Is the node in the neighbour processor?"""
        if not self.is_neighbour(ngh_no):
            return False
        return node_id in self._ngh_ids[str(ngh_no)]

    def get_no_neighbours(self):
        """Get the number of neigbhouring procesor"""
        return len(self._ngh_ids.keys())
        
    def get_no_copies(self, node_id):
        """Get the number of neighbouring processor with a copy of the node"""
        count = 0
        for ngh_no in self._ngh_ids.keys():
            if self.is_neighbour_node(ngh_no, node_id):
                count = count+1
        return count
        
    def get_no_nodes(self, ngh_no):
        """Get the number of overlap nodes in the neighbour processor"""
        if not self.is_neighbour(ngh_no):
            return 0
        return len(self._ngh_ids[str(ngh_no)])

        
    def __iter__(self):
        """Set up an iterator class to loop over the neighbour processor"""
        for key in self._ngh_ids:
            yield self._ngh_ids[key]
            
    def display_neighbour_table(self):
        """Print the neighbour node table
        
        The format is to loop through the table and print each neigbhour number
        and then print all of the overlap nodes in that neighbour
        
        """
        for key in self._ngh_ids:
            print key
            for node_id in self._ngh_ids[key]:
                print "  ", node_id
 