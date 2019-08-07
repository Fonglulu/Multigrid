# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:35:04 2013

@author: stals

This class keeps a record of the local node id and provides some
elementary operations on the id.
"""

class NodeID:
    """ Local node id """
    
    # The local node id is essentially just a counter
    
    _id_no = 0
   
    def __init__(self, id_no=0):
        """Initialise the node id
        
        The id number is set to zero by default
        
        """       
        self._id_no = id_no
        
    def __str__(self):
        """Convert a node id into a string so it can be printed"""
        return  repr(self._id_no)
              
    def __eq__(self, node_id):
        """Check if two node ids are equal"""
        return self._id_no == node_id.get_no()
       
    def __neq__(self, node_id):
        """Check if two node node_ids are not equal"""
        return self._id_no != node_id.get_no()
       
    def __lt__(self, node_id):
        """Check if one node node_id is less than the other"""
        return self._id_no < node_id.get_no()
       
    def __le__(self, node_id):
        """Check if one node id less than or equal to the other"""
        return self._id_no <= node_id.get_no()
       
    def __gt__(self, node_id):
        """Check if one node id greater than the other"""
        return self._id_no > node_id.get_no()
       
    def __ge__(self, node_id):
        """Check if one node id greater than or equal to the other"""
        return self._id_no >= node_id.get_no()       

    def set_no(self, id_no):
        """Set the id number
        
        This should only be used as a help function to other classes
        and it assumed that the number is just some counter
        
        """
        self._id_no = id_no
        
        
    def get_no(self):
        """Get the id number
        
        This should only be used as a help function to other classes
        
        """
        return self._id_no