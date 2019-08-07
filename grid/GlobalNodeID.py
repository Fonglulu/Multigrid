# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 16:35:04 2013

@author: stals

This class keeps a record of the global node id and provides some
elementary operations on the id.

"""

#######################################################################
# mesh_integers_1D
#
# Mesh or combine an integer to form one big integer. This may be used
# to form a global node id. The integer will be of the form x/h and gives
# a unique number to each finite element node.
#
# As this is just a 1D grid, the routine simply returns xh. The routine
# has been included for completeness
#
# Input: integer xh
#
# Output: long integer xh
#
#######################################################################
def mesh_integers_1D(xh):
    """Mesh one integer"""
    return long(xh)
    
#######################################################################
# unmesh_integers_1D
#
# Separate the integer used to form one big integer. 
#
# As this is just a 1D grid, the routine simply returns xh. The routine
# has been included for completeness
#
# Input: integer xh
#
# Output: long integer xh
#
#######################################################################
def unmesh_integers_1D(xh):
    """Unmesh one integer"""
    return xh
    
#######################################################################
# mesh_integers_2D
#
# Mesh or combine two integers to form one big integer. This may be used
# to form a global node id. The integers will be of the form x/h, y/h and 
# give a unique number to each finite element node.
#
# If xh = ....x_21*0^2+x_1*10^1+x_0*10^0 and
#    yh = ....y_2*10^2+y_1*10^1+y_0*10^0 then the meshed integer is
#         ....y_2x_2y_1x_1y_0x_0
#
# Input: positive integer xh
#        positive integer yh
#
# Output: long integer that combines xh and yh
#
#######################################################################
def mesh_integers_2D(xh, yh):
    """Mesh two integers"""
    
    # Make sure the integers are positive 
    assert xh >=0 and yh >= 0, \
            "the input integers should be positive :"+str(xh)+", "+str(yh)

    # Convert the integers into a string
    a = str(xh)
    b = str(yh)
    
    # Pad them out with 0s so both strings are of the same length
    n = max(len(a), len(b))
    Aj = a.rjust(n, '0')
    Bj = b.rjust(n, '0')
    
    # Mesh the two strings together
    D = ''
    for i in range(n):
        D = D + Bj[i]
        D = D + Aj[i]
        
    # Return the result, as a long integer
    return long(D)
    
#######################################################################
# unmesh_integers_2D
#
# Separate the two integers used to form one big integer. 
#
# If xh = ....x_21*0^2+x_1*10^1+x_0*10^0 and
#    yh = ....y_2*10^2+y_1*10^1+y_0*10^0 then the meshed integer is
#         ....y_2x_2y_1x_1y_0x_0
#
# Input: the combined integer
#
# Output: xh, yh
#
#######################################################################
def unmesh_integers_2D(zh):
    """Unmesh two integers"""
    
    # Convert the integer into a string and pad so string length is a
    # multiple of 2
    z = str(zh)
    if len(z)%2 > 0:
        z = z.rjust((len(z)/2+1)*2, '0')
        
    # Easier to step trough the string in reverse order
    zr = z[::-1]
    
    # All of the even entries belong to a and the odd belong to b
    step = len(zr)/2
    ar = ''
    br = ''
    for i in range(step):
        ar = ar + zr[2*i]
        br = br + zr[2*i+1]
        
    # Put the string back in the right order
    A = ar[::-1]
    B = br[::-1]
        
    # Return the result as long integers
    return long(A), long(B)
    


#######################################################################
# mesh_integers_3D
#
# Mesh or combine three integers to form one big integer. This may be used
# to form a global node id. The integers will be of the form x/h, y/h, z/h 
# and give a unique number to each finite element node.
#
# If xh = ....x_21*0^2+x_1*10^1+x_0*10^0,
#    yh = ....y_2*10^2+y_1*10^1+y_0*10^0 and
#    zh = ....z_2*10^2+z_1*10^1+z_0*10^0 then the meshed integer is
#         ....z_2y_2x_2z_1y_1x_1z_0y_0x_0
#
# Input: integer xh
#        integer yh
#        integer zh
#
# Output: long integer that combines xh, y_h and zh
#
#######################################################################
def mesh_integers_3D(xh, yh, zh):
    """Mesh three integers"""
    
    # Make sure the integers are positive 
    assert xh >=0 and yh >= 0 and zh >= 0, \
            "the input integers should be positive"
            
    # Convert the integers into a string
    a = str(xh)
    b = str(yh)
    c = str(zh)
    
    # Pad them out with 0s so all strings are of the same length
    n = max(len(a), len(b), len(c))
    Aj = a.rjust(n, '0')
    Bj = b.rjust(n, '0')
    Cj = c.rjust(n, '0')

    # Mesh the three strings together
    D = ''
    for i in range(n):
        D = D + Cj[i]
        D = D + Bj[i]
        D = D + Aj[i]
        
    # Return the result, as a long integer
    return long(D)

#######################################################################
# unmesh_integers_3D
#
# Separate the three integers used to form one big integer. 
#
# If xh = ....x_21*0^2+x_1*10^1+x_0*10^0,
#    yh = ....y_2*10^2+y_1*10^1+y_0*10^0 and
#    zh = ....z_2*10^2+z_1*10^1+z_0*10^0 then the meshed integer is
#         ....z_2y_2x_2z_1y_1x_1z_0y_0x_0
#
# Input: the combined integer
#
# Output: xh, yh, zh
#
#######################################################################
def unmesh_integers_3D(zh):
    """Unmesh two integers"""
    
    # Convert the integer into a string and pad so string length is a
    # multiple of 3
    z = str(zh)
    if len(z)%3 > 0:
        z = z.rjust((len(z)/3+1)*3, '0')
        
    # Easier to step trough the string in reverse order
    zr = z[::-1]
    
    # All of the even entries below to a and the odd belong to b
    step = len(zr)/3
    ar = ''
    br = ''
    cr = ''
    for i in range(step):
        ar = ar + zr[3*i]
        br = br + zr[3*i+1]
        cr = cr + zr[3*i+2]
        
    # Put the string back in the right order
    A = ar[::-1]
    B = br[::-1]
    C = cr[::-1]
        
    # Return the result
    return long(A), long(B), long(C)
    

###########################################################################

class GlobalNodeID:
    """ Global node id"""
    
    # The global node consists of a large counter and 
    # refinement level
    _id_no = 0
    _level = 0
    
    def __init__(self, id_no = 0, level = 0):
        """Initialise the global node id
        
        The "meshed" number and level are set to 0
        
        """       
        self._id_no = id_no
        self._level = level
        
    def __str__(self):
        """Convert a global node id into a string so it can be printed
        
        The format is the no_level
        
        """
        return  repr(self._id_no)+"_"+repr(self._level)
       
    def __eq__(self, global_id):
        """Check if two global node ids are equal"""
        return self._id_no == global_id.get_no() and \
            self._level == global_id.get_level()
        
    def __neq__(self, global_id):
        """Check if two global node ids are not equal"""
        return self._id_no != global_id.get_no() \
            or self._level != global_id.get_level()
       
    def __lt__(self, global_id):
        """Check if one global node id is less than the other"""
        return (self._id_no < global_id.get_no() \
            and self._level == global_id.get_level()) \
            or self._level < global_id.get_level()
            
    def __le__(self, global_id):
        """Check if one global node id less than or equal to the other"""
        return (self._id_no <= global_id.get_no() \
            and self._level == global_id.get_level()) \
            or self._level < global_id.get_level()
           
    def __gt__(self, global_id):
        """Check if one global node id greater than the other"""
        return (self._id_no > global_id.get_no() \
            and self._level == global_id.get_level()) \
            or self._level > global_id.get_level()
        
    def __ge__(self, global_id):
        """Check if one global node id greater than or equal to the other"""
        return (self._id_no >= global_id.get_no() \
            and self._level == global_id.get_level()) \
            or self._level > global_id.get_level()
   
    def set_no(self, id_no):
        """Set the id number
        
        This should only be used as a help function to other classes
        and it assumed that the number is just some counter
        
        """
        self._id_no = id_no
        
    def set_level(self, level):
        """Set the id level
        
        This should only be used as a help function to other classes
        and it assumed that the number is just some counter
        
        """
        self._level = level
        
    def get_no(self):
        """Get the id number
        
        This should only be used as a help function to other classes
        
        """
        return self._id_no
        
    def get_level(self):
        """Get the id level
        
        This should only be used as a help function to other classes
        
        """
        return self._level
