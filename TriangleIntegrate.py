# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:35:04 2013

@author: stals

This class offers a number of numerical quadrature routines to integrate
a function over a triangle.

Currently only offers quadrature routines that are exact for linear 
functions

There are generally two types of rountines, one that integrates a 
polynomial*polynomial and one that integrates a function*polynomial.
"""


#######################################################################
# calc_area
#
# Calculate the area of the triangle (coord1, coord2, coord3).
#
# It is assumed that the length of all of the coordinates is greather 
# than 1.
#
# Input: array coord1
#        array coord2
#        array coord3
#
# Output: float area
#
#######################################################################
def calc_area(coord1, coord2, coord3):
    """Calculate the area of the triangle"""

    from math import fabs
    return 0.5 * fabs((coord2[0]-coord1[0])*(coord3[1]-coord1[1])
        - (coord3[0]-coord1[0])*(coord2[1]-coord1[1]))     

######################################################################
# linear_integrate
#
# Estimate the integral of uv, where u and v are polynomials, using 
# a quadrature rule which is exact for polynomials of degree 1.
#
#  See Cleas Johnson: Numerical solution of partial differential 
# equations by the finite element method, Cambridge, 1987
#
# It is assumed that node1, node2 and node3 are all two dimensional 
# nodes.
#
# Input: A two dimensional polynomial poly1
#        A two dimensional polynomial poly2
#        Node node1
#        Node node2
#        Node node3
#
# Output: float approximate integral
#
#######################################################################
def linear_integrate(poly1, poly2, node1, node2, node3): 
    """Evaluate the approximate integral over the given triangle"""

    # Import the appropriate information from other modules
    from numpy import array
    
    # Check the coordinates are two dimensional
    assert node1.get_dim() == 2 and node2.get_dim() == 2 \
        and node3.get_dim() == 2, \
            "the triangle coordinates must be two dimensional"
            
    # Get the coordinates of the 3 vertices
    coord1 = array(node1.get_coord())
    coord2 = array(node2.get_coord())
    coord3 = array(node3.get_coord())
    
    # Find the point at the centre of the triangle
    centre_point = (coord1+coord2+coord3)/3.0
    x = centre_point[0]
    y = centre_point[1]
    
    # Apply the quadrature rule    
    return  poly1.eval(x, y) * poly2.eval(x, y) * \
        calc_area(coord1, coord2, coord3)
        

######################################################################
# linear_func_integrate
#  Estimate the integral of fv, where f is a function and v is a 
# polynomials, using  a quadrature rule which is exact for 
# polynomials of degree 1.
#
#  See Cleas Johnson: Numerical solution of partial differential 
# equations by the finite element method, Cambridge, 1987
#
# It is assumed that node1, node2 and node3 are all two dimensional 
# nodes.
#
# Input: A two dimensional polynomial poly1
#        A two dimensional polynomial poly2
#        Node node1
#        Node node2
#        Node node3
#
# Output: float approximate integral
#
#######################################################################
def linear_func_integrate(rhs, poly, node1, node2, node3): 
    """Evaluate the approximate integral over the given triangle"""
    
    # Import the appropriate information from other modules
    from numpy import array
    
    # Check the coordinates are two dimensional
    assert node1.get_dim() == 2 and node2.get_dim() == 2 \
     and node3.get_dim() == 2, \
            "the triangle coordinates must be two dimensional"

    # Get the coordinates of the 3 vertices
    coord1 = array(node1.get_coord())
    coord2 = array(node2.get_coord())
    coord3 = array(node3.get_coord())
    
    # Find the point at the centre of the triangle
    centre_point = (coord1+coord2+coord3)/3.0
    x = centre_point[0]
    y = centre_point[1]
    
    # Apply the quadrature rule    
    return  rhs(centre_point) * poly.eval(x, y) * \
        calc_area(coord1, coord2, coord3)