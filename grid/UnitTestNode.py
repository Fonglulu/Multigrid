# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:35:04 2013

@author: stals

This is a collection of unit tests for the node classes

It has only been given superficial documentation
"""

# import the appropriate information from other classes

import unittest
from GlobalNodeID import GlobalNodeID
from NodeID import NodeID
from Node import Node
from NodeTable import NodeTable

class TestNodeID(unittest.TestCase):
    """ Set up a suit of tests for the NodeID class"""
   
    
    def set_up(self):
        id1 = NodeID()
        id2 = NodeID()
        id1.set_no(1)
        id2.set_no(0)
        return [id1, id2]
        
    def test_set_get_no(self):
        id1, id2 = self.set_up()
        self.assertEqual(id1.get_no(), 1)
        
    def test_eq(self):
        id1, id2 = self.set_up()
        self.assertEqual(id1 == id2, False)

    def test_neq(self):
        id1, id2 = self.set_up()
        self.assertEqual(id1 != id2, True)
        
    def test_lt(self):
        id1, id2 = self.set_up()
        self.assertEqual(id1 < id2, False)
        
    def test_le(self):
        id1, id2 = self.set_up()
        self.assertEqual(id1 <= id2, False)

    def test_gt(self):
        id1, id2 = self.set_up()
        self.assertEqual(id1 > id2, True)

    def test_ge(self):
        id1, id2 = self.set_up()
        self.assertEqual(id1 >= id2, True)
        
        
class TestGlobalID(unittest.TestCase):
    """ Set up a suit of tests for the GlobalNodeID class"""
    
    def set_up(self):
        global_id1 = GlobalNodeID()
        global_id2 = GlobalNodeID()
        global_id1.set_no(1)
        global_id1.set_level(1)
        global_id2.set_no(0)
        global_id2.set_level(1)
        return global_id1, global_id2
        
    def test_set_get_no(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1.get_no(), 1)

    def test_set_get_level(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1.get_level(), 1)
        
    def test_eq(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1 == global_id2, False)

    def test_neq(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1 != global_id2, True)
        
    def test_lt(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1 < global_id2, False)
        
    def test_le(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1 <= global_id2, False)

    def test_gt(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1 > global_id2, True)

    def test_ge(self):
        global_id1, global_id2 = self.set_up()
        self.assertEqual(global_id1 >= global_id2, True)



        
class TestNode(unittest.TestCase):
    """ Set up a suit of tests for the Node class"""
   
    def set_up(self):
        node1 = Node()
        node2 = Node()
        node_id = NodeID()
        global_id = GlobalNodeID()
        node_id.set_no(2)
        coord = [1, 4]
        global_id.set_no(1)
        global_id.set_level(2)
        node1.set_node_id(node_id)
        node1.set_global_id(global_id)
        node1.set_value(-1.5)
        node1.set_slave(False)
        node1.set_load(0.1)
        node1.set_coord(coord)
        node_id.set_no(3)
        global_id.set_no(2)
        global_id.set_level(2)
        coord = [-5 , 7]
        node2.set_node_id(node_id)
        node2.set_global_id(global_id)
        node2.set_value(-2.5)
        node2.set_slave(True)
        node2.set_load(0.5)
        node2.set_coord(coord)
 
        return [node1, node2]
        
    def test_set_get_id(self):
        node1, node2 = self.set_up()
        node_id = NodeID()
        node_id.set_no(2)
        self.assertEqual(node1.get_node_id(), node_id)
        
    def test_set_get_global_id(self):
        node1, node2 = self.set_up()
        global_id = GlobalNodeID()
        global_id.set_no(1)
        global_id.set_level(2)
        self.assertEqual(node1.get_global_id(), global_id)
        
    def test_set_get_coord(self):
        node1, node2 = self.set_up()
        coord = [1, 4]
        self.assertEqual(node1.get_coord(), coord)
        
    def test_set_get_load(self):
        node1, node2 = self.set_up()
        self.assertEqual(node1.get_load(), 0.1)

    def test_set_get_value(self):
        node1, node2 = self.set_up()
        self.assertEqual(node1.get_value(), -1.5)

    def test_set_get_slave(self):
        node1, node2 = self.set_up()
        self.assertEqual(node1.get_slave(), False)
 
        
class TestNodeTable(unittest.TestCase):
    """ Set up a suit of tests for the NodeTable class"""
   
    def set_up(self):
        node_table = NodeTable()
        node1 = Node()
        node2 = Node()
        node_id = NodeID()
        global_id = GlobalNodeID()
        node_id.set_no(2)
        coord = [1, 4]
        global_id.set_no(1)
        global_id.set_level(2)
        node1.set_node_id(node_id)
        node1.set_global_id(global_id)
        node1.set_value(-1.5)
        node1.set_slave(False)
        node1.set_load(0.1)
        node1.set_coord(coord)
        node_id.set_no(3)
        global_id.set_no(2)
        global_id.set_level(2)
        coord = [-5 , 7]
        node2.set_node_id(node_id)
        node2.set_global_id(global_id)
        node2.set_value(-2.5)
        node2.set_slave(True)
        node2.set_load(0.5)
        node2.set_coord(coord)
        node_table.add_node(node1)
        node_table.add_node(node2)
 
        return node_table
        
    def test_set_get_id(self):
        node_table = self.set_up()
        global_id = GlobalNodeID()
        true_id = NodeID()
        global_id.set_no(1)
        global_id.set_level(2)
        true_id.set_no(2)
        node_id = node_table.get_local_id(global_id)
        self.assertEqual(node_id, true_id)
        
    def test_set_get_global_id(self):
        node_table = self.set_up()
        true_global_id = GlobalNodeID()
        node_id = NodeID()
        true_global_id.set_no(1)
        true_global_id.set_level(2)
        node_id.set_no(2)
        global_id = node_table.get_global_id(node_id)
        self.assertEqual(global_id, true_global_id)
        
    def test_set_get_coord(self):
        node_table = self.set_up()
        node_id = NodeID()
        node_id.set_no(2)
        coord = [1, 4]
        self.assertEqual(node_table.get_coord(node_id), coord)
        
    def test_set_get_load(self):
        node_table = self.set_up()
        node_id = NodeID()
        node_id.set_no(2)
        self.assertEqual(node_table.get_load(node_id), 0.1)

    def test_set_get_value(self):
        node_table = self.set_up()
        node_id = NodeID()
        node_id.set_no(2)
        self.assertEqual(node_table.get_value(node_id), -1.5)

    def test_set_get_slave(self):
        node_table = self.set_up()
        node_id = NodeID()
        node_id.set_no(2)
        self.assertEqual(node_table.get_slave(node_id), False)
        
    def test_set_get_dim(self):
        node_table = self.set_up()
        self.assertEqual(node_table.get_dim(), 2) 
        
    def test_is_in(self):
        node_table = self.set_up()
        node_id = NodeID()
        node_id.set_no(1)
        self.assertEqual(node_table.is_in(node_id), False)
        
    def test_delete(self):
        node_table = self.set_up()
        node_id = NodeID()
        node_id.set_no(2)
        node_table.delete_node(node_id)
        self.assertEqual(node_table.is_in(node_id), False)
        
    def test_count(self):
        node_table = self.set_up()
        self.assertEqual(node_table.get_no_nodes(), 2)

if __name__ == '__main__':
    unittest.main()
