# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:35:04 2013

@author: stals

This is a collection of unit tests for the neighbour node class

It has only been given superficial documentation
"""

# Import the appropriate information from other classes
import unittest
from NghNodes import NghNodes
from NodeID import NodeID

class TestNghNode(unittest.TestCase):
    """ Set up a suit of tests for the NodeID class"""
   
    
    def set_up(self):
        id1 = NodeID(1)
        id2 = NodeID(2)
        id3 = NodeID(3)
        ngh_nodes = NghNodes()
        ngh_nodes.add_neighbour_node(1, id1)
        ngh_nodes.add_neighbour_node(1, id2)
        ngh_nodes.add_neighbour_node(2, id3)
        ngh_nodes.add_neighbour_node(2, id2)
        return ngh_nodes
        
    def test_add_neighbour(self):
        ngh_nodes = self.set_up()
        id2 = NodeID(2)
        self.assertTrue(ngh_nodes.get_no_neighbours() == 2 \
                and ngh_nodes.get_no_copies(id2) == 2 \
                and ngh_nodes.get_no_nodes(1) == 2)
                
    def test_delete_neighbour(self):
        ngh_nodes = self.set_up()
        id2 = NodeID(2)
        ngh_nodes.delete_neighbour_node(1, id2)
        self.assertTrue(ngh_nodes.get_no_neighbours() == 2 \
                and ngh_nodes.get_no_copies(id2) == 1 \
                and ngh_nodes.get_no_nodes(1) == 1)
 

if __name__ == '__main__':
    unittest.main()
