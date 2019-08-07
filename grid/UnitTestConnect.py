# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:32:19 2013

@author: stals

This is a collection of unit tests for the connect classes

It has only been given superficial documentation
"""

# import the appropriate information from other classes
import unittest
from copy import copy
from NodeID import NodeID
from ConnectStar import ConnectStar, ConnectEnd
from ConnectTable import ConnectTable


        
class TestConnectEnd(unittest.TestCase):
    """ Set up a suit of tests for the ConnectEnd class"""
   
    def set_up(self):
        node_id1 = NodeID(1)
        connect_end = ConnectEnd()
        connect_end.set_id(node_id1)
        connect_end.set_value(1.0)
        return connect_end
        
    def test_set_get_id(self):
        connect_end = self.set_up()
        node_id = NodeID(1)
        self.assertEqual(connect_end.get_id(), node_id)
        
    def test_set_get_value(self):
        connect_end = self.set_up()
        self.assertEqual(connect_end.get_value(), 1.0)
        
        
class TestconnectStar(unittest.TestCase):
    """ Set up a suit of tests for the connectStar class"""
   
    def set_up(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        node_id3 = NodeID(3)
        connect_star = ConnectStar(node_id1)
        connect_star.set_id(node_id1)
        connect_star.add_connection(node_id2)
        connect_star.set_value(node_id2, 2.0)
        connect_star.add_connection(node_id3)
        connect_star.set_value(node_id3, 3.0)
        return connect_star
        
    def test_set_get_id(self):
        connect_star = self.set_up()
        node_id = NodeID(1)
        self.assertEqual(connect_star.get_id(), node_id)
        
        
    def test_set_get_value(self):
        node_id2 = NodeID(2)
        connect_star = self.set_up()
        self.assertEqual(connect_star.get_value(node_id2), 2.0)
        
    def test_get_no_endpoints(self):
        connect_star = self.set_up()
        self.assertEqual(connect_star.get_no_endpoints(), 2)
        
    def test_delete(self):
        connect_star = self.set_up()
        node_id2 = NodeID(2)
        connect_star.delete_connection(node_id2)
        self.assertEqual(connect_star.get_no_endpoints(), 1)
        
        
        
class TestconnectTable(unittest.TestCase):
    """ Set up a suit of tests for the connectTable class"""
   
    def set_up(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        node_id3 = NodeID(3)
        connect_table = ConnectTable()
        connect_table.add_connection(node_id1, node_id2, 1.0)
        connect_table.add_connection(node_id1, node_id3, 2.0)
        connect_table.add_connection(node_id2, node_id3, 3.0)
        return connect_table
        
    def test_set_get_value(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        connect_table = self.set_up()
        self.assertEqual(connect_table.get_matrix_value(node_id1, node_id2), \
            1.0)
        
    def test_get_no_endpoints(self):
        node_id1 = NodeID(1)
        connect_table = self.set_up()
        self.assertEqual(connect_table.get_no_endpoints(node_id1), 2)
        
    def test_delete(self):
        connect_table = self.set_up()
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        connect_table.delete_connection(node_id1, node_id2)
        self.assertEqual(connect_table.get_no_connections(), 2)
        
if __name__ == '__main__':
    unittest.main()
 