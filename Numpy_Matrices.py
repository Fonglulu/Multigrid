#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 09:07:36 2019

@author: shilu
"""

import numpy as np
from copy import copy
from functions import Exy

def Numpy_L(n):
    
    pass



class Triangle_container:
    """ Triagles asscociate with data points"""
    
    _data = np.zeros((1,3))
    
    _coord = np.array([])
    
    _vert_id  = []
    
    
    def __init__(self, data = np.zeros((1,3)), coord =np.array([])):
        
        
        self._data = data
        self._coord = coord
        #self._veverticesrtices = 
        
    def add_data(self, a_data):
        
        self._data = np.concatenate((self._data, a_data),axis=0)
        
        #return self._data
    
    def set_coord(self, vert):
        #print 1
        ''' Set up the coordinates of the three vertices'''
        self._coord = vert[:]
    
   # def get_triangle(coord):
        
        #return 
    
#    def set_vertex(self, index1, index2, index3):
#        
#        ''' Set up the index of the three vertices'''
#        
#        self._vertices = [index1, index2, index3]
#        

        
        
        
    
    
class Triangle_table:
    
    
    _triangles = np.array([])
    
    _vertices = []
    
    def __init__(self, triangles = np.array([]), vertices = []):
        
        self._triangles = triangles
        
        self._vertices = vertices
        
    def add_vertices(self, coord1, coord2, coord3):
        #print '1' 
        
        vert = np.array([coord1, coord2, coord3])
        
        self._vertices.append(vert)
        
        
    def get_triangle(self, coord1, coord2, coord3):
        
        
        vert = np.array([coord1, coord2, coord3])
        
        #print self._triangles[0],'1'
        
        for triangle in self._triangles:
            
             #print triangle
             
             if np.all(triangle._coord == vert):
        
        
                 return triangle #for triangle in self._triangles if np.all(triangle._coord == vert)]
    
    
        
    def create_triangle(self,  coord1, coord2, coord3):
           
            triangle = Triangle_container()
            #print '1'
            #print triangle._data
            #triangle.add_data(data)
            
            vert = np.array([coord1, coord2, coord3])
            
            triangle.set_coord(vert)
            
            return triangle
            
            
            
        
    def add_triangle(self, triangle):
            
            self._triangles = np.append(self._triangles, copy(triangle))
            

        
        
        
            
            
         
            
def pre_processing(i):
    

    grid = Triangle_table()
    
    n = 2**i+1
    
    h =1/float(n-1)
    
    
    x = np.linspace(0, 1.0,6)
    y = np.linspace(0, 1.0,6)
    X, Y = np.meshgrid(x,y)
    
    model_func =Exy(X,Y)
    
    
    value = model_func.flatten()
    coordx = X.flatten()
    coordy = Y.flatten()
    Coord = zip(coordx, coordy)
    
    for i in range(len(Coord)):
        print i
 
        data = Coord[i]
        
    
        if (data[0] !=0) and (data[0]!=1) and (data[1]!= 0 ) and (data[1] != 1): 
            #print data, 'data'
            num_x = int(data[0]/float(h))
    
            num_y = int(data[1]/float(h))
        
        
            left_upper = (num_x*h, (num_y+1)*h)
            
            lu_id = num_x + n*(num_y+1)
            
            print left_upper, lu_id
       
            left_down = (num_x*h, num_y*h)
        
            right_upper = ((num_x+1)*h, (num_y+1)*h)
        
            right_down = ((num_x+1)*h, num_y*h)
            #print data, left_down, right_down, right_upper
            
            # If the data is in the lower triangle
            if In_triangle(left_down, right_down, right_upper, data):
            
                #print left_down, right_down, right_upper, 'lower'
            
           
            
                vert = np.array([left_down, right_down, right_upper])
                #print vert 
            
                # If the group of three vertices has been not been stored 
                if filter(lambda x: np.all(x == vert), grid._vertices) ==[]:
                    #print '2'
                    grid.add_vertices(left_down, right_down, right_upper)
                
                    tri = grid.create_triangle( left_down, right_down, right_upper)
                
                    #grid.add_triangle(tri)
                    
                    data_value = value[i]
                   
                    data_3D = data+(data_value,)
                    
                    data_3D = np.array([data_3D])
                    #print tri._data, data_3D
                    tri.add_data(data_3D)
                    
                    #print tri._coord, tri._data
                    
                    grid.add_triangle(tri)
                    #print grid._triangles[0]._data, 'data'
                    
                else:
                
                        tri = grid.get_triangle(left_down, right_down, right_upper)
                
                        #print '3'
                        
                        data_value = value[i]
                    
                        data_3D = data+(data_value,)
                    
                        data_3D = np.array([data_3D])
                        
                        tri.add_data(data_3D)
                    
                
                
                
                
                
                
                
                
                
                
            # If the data is in the upper triangle 
            else:
                
                vert = np.array([left_down, left_upper, right_upper])
                #print vert , 'upper'
                if filter(lambda x: np.all(x == vert), grid._vertices) ==[]:
            
                    tri = grid.create_triangle(left_down, left_upper, right_upper)
            
                    
                    
                    
                    data_value = value[i]
                    
                    data_3D = data+(data_value,)
                    
                    data_3D = np.array([data_3D])
                    
                    tri.add_data(data_3D)
                    
                    grid.add_triangle(tri)
                

                
                else:
                
                    tri = grid.get_triangle(left_down, right_upper, right_upper)
                    
                    data_value = value[i]
                    
                    data_3D = data+(data_value,)
                    
                    data_3D = np.array([data_3D])
                    
                    tri.add_data(data_3D)

        
    return grid
            
        
        
        
        
    



def Polynomial_eval(node1, node2, node3, data_coord):
    
    from math import fabs
    #print data_coord
    x1 = node1[0]
    y1 = node1[1]
    x2 = node2[0]
    y2 = node2[1]
    x3 = node3[0]
    y3 = node3[1]
    
    division = (y1-y2)*(x2-x3)-(y2-y3)*(x1-x2);
    #print division
    
    
    assert fabs(division) > 1.0E-12, "divide by zero"
    
    const = (x3*y2 - y3*x2)/float(division)
    xcoe = (y3-y2)/float(division)
    ycoe = (x2-x3)/float(division)
    
    return data_coord[0]*xcoe+data_coord[1]*ycoe+const


def In_triangle(node1, node2, node3, data_coord):
    
    value1 = Polynomial_eval(node1, node2, node3, data_coord)
    value2 = Polynomial_eval(node2, node1, node3, data_coord)
    value3 = Polynomial_eval(node3, node2, node1, data_coord)
    
    return  (value1 >=0.0 and value2>=0.0 and value3 >=0.0)

    
            
        
        
        
        
        
        
