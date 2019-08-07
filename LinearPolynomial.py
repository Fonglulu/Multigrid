# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:35:04 2013

@author: stals

This class keeps record of a linear polynomial in two dimensions.
It currently only offers minimal functionality, additional polynomial
opertations, such as addition and scalar multiplication, could also 
be included.
"""

        
#######################################################################

class LinearPolynomial:
    """ A linear polynomial in two dimensions"""
    
    # constant coefficient
    _const = 0.0
    
    # x coefficient
    _x = 0.0
    
    # y coefficient
    _y = 0.0
   
   
    def __init__(self, const = 0.0, x = 0.0, y = 0.0):
        """Initialise a polynomial
        
        By default, the polynomial is the zero polynomial
        """
        self._const = const
        self._x = x
        self._y = y

        
    def __str__(self):
        """Convert a polynomial into a string so it can be printed
        
        The format is constant coefficient + x coefficient * x 
        + y coefficient * y
        
        """
        return str(self._const) + " + " + str(self._x) + "x + " \
        + str(self._y) + "y"
    
     
    def set_const(self, const):
        """Set constant coefficient"""
        self._const = const
    
    def set_x(self, x):
        """Set the x coefficient"""
        self._x = x
        
    def set_y(self, y):
        """Set the y coefficient"""
        self._y = y

    def get_const(self):
        """Get the constant coefficient"""
        return self._const
        
    def get_x(self):
        """Get the x coefficient"""
        return self._x
    
    def get_y(self):
        """Get the y coefficient"""
        return self._y
        
    def eval(self, x, y): 
        """evaluate the polynomial at the posn (x, y)"""
        return self._const + self._x * x + self._y * y
        
    def dx(self): 
        """differentiate the polynomial in the x direction"""
        poly = LinearPolynomial(self._x, 0.0, 0.0)
        return poly
        
    def dy(self): 
        """differentiate the polynomial in the y direction"""
        poly = LinearPolynomial(self._y, 0.0, 0.0)
        return poly
        