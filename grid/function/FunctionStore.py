# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 08:44:30 2013

@author: stals

This modules lists all of the functions used to define the example problems.
That is, it lists the boundary conditions and right hand side functions.

Any new functions added to the list should also be included in the FunctionStore
class. The purpose of the class is to get around the issue of not being able
to communicate a function through MPI. As each process will have access to the
same list of functions given in FunctionStore, MPI can (in essence) 
communicate the position of the function in the store. By knowing the position,
each process can then determine what function is being communicated.

Note that with Python this trick (or hack) is not strictly needed, Python's
pickle command can be used to communicate a function. Nevertheless, I have
decided to continue to use the FunctionStore as pickle is not, generally,
available in other languages and the use of pickle is inefficient. 

"""

############################################################################
# We firstly define some functions that will be used in the test problems.
############################################################################

# Import appropriate information from other modules
from numpy import pi, sin, exp, cos

def zero(x):
    """Define the function 0"""
    return 0

def plain(x):
    
    return x[0]-x[0]+1.0

def l2(x):
    
    return 2*x[0] #x[0]**2+x[1]**2 

def l3(x):
    
    return  (x[0]+1)**3+x[1]**3

def linear(x):
    
    return x[0]+x[1]

def linear_x(x):
    
    return 2*x[0]+x[1]
    
def exp_soln(x):
    """Define the function exp(x_0)exp(x_1)"""
    return exp(x[0])*exp(x[1])
    
def exp_rhs(x):
    """Define the function -2*exp(x_0)exp(x_1)"""
    return -2.0*exp(x[0])*exp(x[1])
    
def sin_soln(x):
    """Define the function sin(pi x_0)sin(pi x_1)"""
    return sin(pi*x[0])*sin(pi*x[1])
    
def sin_rhs(x):
    """Define the function 2pi^2 sin(pi x_0)sin(pi x_1)"""
    return 2.0*pi*pi*sin(pi*x[0])*sin(pi*x[1])
    
def sqr_soln(x):
    """Define the function x_0^2+x_1^2"""
    return x[0]*x[0]+x[1]*x[1]
    
def sqr_rhs(x):
    """Define the function -4"""
    return -4.0

def sin4(x):
    
    return sin(4*pi*x[0])*sin(4*pi*x[1])

def exy(x):
    return exp(3/((x[0]-0.5)**2+(x[1]-0.5)**2+1))


def cos2(x):
    
    return cos(2*pi*x[0])*cos(2*pi*x[1])

############################################################################
# We now include the functions in a list that will be available on all
# processors
############################################################################
class FunctionStore:
    """Store the functions used in the test problems"""
    
    # Create an empty list
    _func_list = list()
    
    def __init__(self):
        """Initialise the function store
        
        The functions used in the test problems should be appended here"""
        self._func_list = [zero]
        self._func_list.append(plain)
        self._func_list.append(l2)
        self._func_list.append(l3)
        self._func_list.append(linear)
        self._func_list.append(exp_soln)
        self._func_list.append(exp_rhs)
        self._func_list.append(sin_soln)
        self._func_list.append(sin_rhs)
        self._func_list.append(sqr_soln)
        self._func_list.append(sqr_rhs)
        self._func_list.append(sin4)
        self._func_list.append(exy)
        self._func_list.append(cos2)
        
    def is_in_store(self, func):
        """Is a given function in the store?"""
        return func in self._func_list
        
    def get_index(self, func):
        """Get the position of the function in the store"""
        assert self.is_in_store(func), \
            "function not in store "
        return self._func_list.index(func)
        
    def get_function(self, index):
        """Get the function that is in the given position in the store"""
        assert index >= 0 and index < len(self._func_list), \
            "index out of range "
        return self._func_list[index]
    
    
