# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:32:19 2013

@author: stals

This is a collection of unit tests for the edge classes

It has only been given superficial documentation
"""

# import the appropriate information from other classes
import unittest

from copy import copy
from NodeID import NodeID
from Edge import Edge, DomainSet, RefineSet
from function.FunctionStore import zero
from EdgeStar import EdgeStar, EdgeEnd
from EdgeTable import EdgeTable

class TestEdge(unittest.TestCase):
    """ Set up a suit of tests for the Edge class"""
   
    def set_up(self):
        endpt1 = NodeID(1)
        endpt2 = NodeID(2)
        edge = Edge()
        edge.set_endpt1(endpt1)
        edge.set_endpt2(endpt2)
        edge.set_location(DomainSet.boundary)
        edge.set_boundary_function(zero)
        edge.set_refine_type(RefineSet.base_edge)
        edge.set_error_indicator(2.0)
        edge.set_refine_level(1)
        return edge
        
    def test_set_get_endpt1(self):
        edge = self.set_up()
        node_id = NodeID(1)
        self.assertEqual(edge.get_endpt1(), node_id)
        
    def test_set_get_endpt2(self):
        edge = self.set_up()
        node_id = NodeID(2)
        self.assertEqual(edge.get_endpt2(), node_id)
        
    def test_set_get_boundary_location(self):
        edge = self.set_up()
        self.assertEqual(edge.get_location(), DomainSet.boundary)
        
    def test_set_get_boundary_function(self):
        edge = self.set_up()
        self.assertEqual(edge.get_boundary_function(), zero)
        
    def test_set_get_error_indicator(self):
        edge = self.set_up()
        self.assertEqual(edge.get_error_indicator(), 2.0)
        
    def test_set_get_refine_level(self):
        edge = self.set_up()
        self.assertEqual(edge.get_refine_level(), 1)
        
    def test_set_get_refine_type(self):
        edge = self.set_up()
        self.assertEqual(edge.get_refine_type(), RefineSet.base_edge)
        
class TestEdgeEnd(unittest.TestCase):
    """ Set up a suit of tests for the EdgeEnd class"""
   
    def set_up(self):
        node_id1 = NodeID(1)
        edge_end = EdgeEnd()
        edge_end.set_id(node_id1)
        edge_end.set_location(DomainSet.boundary)
        edge_end.set_boundary_function(zero)
        edge_end.set_refine_type(RefineSet.base_edge)
        edge_end.set_error_indicator(2.0)
        edge_end.set_refine_level(1)
        return edge_end
        
    def test_set_get_id(self):
        edge_end = self.set_up()
        node_id = NodeID(1)
        self.assertEqual(edge_end.get_id(), node_id)
        
    def test_set_get_boundary_location(self):
        edge_end = self.set_up()
        self.assertEqual(edge_end.get_location(), DomainSet.boundary)
        
    def test_set_get_boundary_function(self):
        edge_end = self.set_up()
        self.assertEqual(edge_end.get_boundary_function(), zero)
        
    def test_set_get_error_indicator(self):
        edge_end = self.set_up()
        self.assertEqual(edge_end.get_error_indicator(), 2.0)
        
    def test_set_get_refine_level(self):
        edge_end = self.set_up()
        self.assertEqual(edge_end.get_refine_level(), 1)
        
    def test_set_get_refine_type(self):
        edge_end = self.set_up()
        self.assertEqual(edge_end.get_refine_type(), RefineSet.base_edge)
        
class TestEdgeStar(unittest.TestCase):
    """ Set up a suit of tests for the EdgeStar class"""
   
    def set_up(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        node_id3 = NodeID(3)
        edge_star = EdgeStar()
        edge_star.set_id(node_id1)
        edge_star.add_endpoint(node_id2)
        edge_star.set_location(node_id2, DomainSet.boundary)
        edge_star.set_boundary_function(node_id2, zero)
        edge_star.set_refine_type(node_id2, RefineSet.base_edge)
        edge_star.set_error_indicator(node_id2, 2.0)
        edge_star.set_refine_level(node_id2, 1)
        edge_star.add_endpoint(node_id3)
        edge_star.set_location(node_id3, DomainSet.interior)
        edge_star.set_boundary_function(node_id3, zero)
        edge_star.set_refine_type(node_id3, RefineSet.not_base_edge)
        edge_star.set_error_indicator(node_id3, -10.0)
        edge_star.set_refine_level(node_id3, 5)
        return edge_star
        
    def test_set_get_id(self):
        edge_star = self.set_up()
        node_id = NodeID(1)
        self.assertEqual(edge_star.get_id(), node_id)
        
    def test_set_get_boundary_location(self):
        node_id2 = NodeID(2)
        edge_star = self.set_up()
        self.assertEqual(edge_star.get_location(node_id2), DomainSet.boundary)
        
    def test_set_get_boundary_function(self):
        node_id2 = NodeID(2)
        edge_star = self.set_up()
        self.assertEqual(edge_star.get_boundary_function(node_id2), zero)
        
    def test_set_get_error_indicator(self):
        node_id2 = NodeID(2)
        edge_star = self.set_up()
        self.assertEqual(edge_star.get_error_indicator(node_id2), 2.0)
        
    def test_set_get_refine_level(self):
        node_id2 = NodeID(2)
        edge_star = self.set_up()
        self.assertEqual(edge_star.get_refine_level(node_id2), 1)
        
    def test_set_get_refine_type(self):
        node_id2 = NodeID(2)
        edge_star = self.set_up()
        self.assertEqual(edge_star.get_refine_type(node_id2), \
            RefineSet.base_edge)
        
    def test_get_no_endpoints(self):
        edge_star = self.set_up()
        self.assertEqual(edge_star.get_no_endpoints(), 2)
        
    def test_delete(self):
        edge_star = self.set_up()
        node_id2 = NodeID(2)
        edge_star.delete_endpoint(node_id2)
        self.assertEqual(edge_star.get_no_endpoints(), 1)
        
        
        
class TestEdgeTable(unittest.TestCase):
    """ Set up a suit of tests for the EdgeTable class"""
   
    def set_up(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        node_id3 = NodeID(3)
        edge1 = Edge(node_id1, node_id2)
        edge1.set_location(DomainSet.boundary)
        edge1.set_boundary_function(zero)
        edge1.set_refine_type(RefineSet.base_edge)
        edge1.set_error_indicator(2.0)
        edge1.set_refine_level(1)
        edge2 = Edge(node_id1, node_id3)
        edge2.set_location(DomainSet.interior)
        edge2.set_boundary_function(zero)
        edge2.set_refine_type(RefineSet.not_base_edge)
        edge2.set_error_indicator(-10.0)
        edge2.set_refine_level(5)
        edge3 = copy(edge2)
        edge3.set_endpt1(node_id2)
        edge_table = EdgeTable()
        edge_table.add_edge_all(edge1)
        edge_table.add_edge_all(edge2)
        edge_table.add_edge_all(edge3)
        return edge_table
        
    def test_set_get_edge(self):
        edge_table = self.set_up()
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        edge = edge_table.get_edge(node_id1, node_id2)
        self.assertTrue(edge.get_endpt1() == node_id1 \
            and edge.get_endpt2() == node_id2\
            and edge.get_location() == DomainSet.boundary\
            and edge.get_boundary_function() == zero \
            and edge.get_refine_type() == RefineSet.base_edge\
            and edge.get_error_indicator() == 2.0\
            and edge.get_refine_level() == 1)

    def test_set_get_boundary_location(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        edge_table = self.set_up()
        self.assertEqual(edge_table.get_location(node_id1, node_id2), \
            DomainSet.boundary)
        
    def test_set_get_boundary_function(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        edge_table = self.set_up()
        self.assertEqual(edge_table.get_boundary_function(node_id1, node_id2),\
            zero)
        
    def test_set_get_error_indicator(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        edge_table = self.set_up()
        self.assertEqual(edge_table.get_error_indicator(node_id1, node_id2), \
            2.0)
        
    def test_set_get_refine_level(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        edge_table = self.set_up()
        self.assertEqual(edge_table.get_refine_level(node_id1, node_id2), 1)
        
    def test_set_get_refine_type(self):
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        edge_table = self.set_up()
        self.assertEqual(edge_table.get_refine_type(node_id1, node_id2), \
            RefineSet.base_edge)
        
    def test_get_no_endpoints(self):
        node_id1 = NodeID(1)
        edge_table = self.set_up()
        self.assertEqual(edge_table.get_no_endpoints(node_id1), 2)
        
    def test_delete(self):
        edge_table = self.set_up()
        node_id1 = NodeID(1)
        node_id2 = NodeID(2)
        edge_table.delete_edge(node_id1, node_id2)
        self.assertEqual(edge_table.get_no_edges(), 2)


if __name__ == '__main__':
    unittest.main()
 
