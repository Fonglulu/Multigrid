#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 09:07:36 2019

@author: shilu
"""

import numpy as np
from copy import copy
from functions import Exy, Xexy, Yexy, XYexy
from copy import deepcopy
from scipy.sparse import  lil_matrix
import time 
def Numpy_L(n):
    
    pass



class Triangle_container:
    """ Triagles asscociate with data points"""
    
    _data = np.zeros((1,3))
    
    _coord = np.array([])
    
    _node_id = set()
    
    _int_node_id = set()

    def __init__(self, data = np.zeros((1,3)), coord =np.array([]), node_id=set(), int_node_id=set()):
        
        
        
#        print self._data
#        print self._coord
#        print self._node_id
        self._data = data
        self._coord = coord
        self._node_id = deepcopy(node_id)
        self._int_node_id =  deepcopy(int_node_id)
        
        
    def add_data(self, a_data):

        if np.all(self._data[0] ==np.zeros((3,))):
            
            self._data = np.delete(self._data,0 ,0)
            
        self._data = np.concatenate((self._data, a_data),axis=0)
        
        
        
        #return self._data
    
    def set_coord(self, vert):

        ''' Set up the coordinates of the three vertices'''
        self._coord = vert[:]
        
    def add_node_id(self, node_id):

        self._node_id.add(node_id)
        
        
    def add_int_node_id(self, int_node_id):
        
        self._int_node_id.add(int_node_id)

    


    
    
class Triangle_table:
    
    
    _triangles = {}
    
    _vertices = []
    
    def __init__(self, triangles = {}, vertices = []):
        
        self._triangles = deepcopy(triangles)
        
        self._vertices = vertices
        
    def add_vertices(self, coord1, coord2, coord3):
        #print '1' 
        
        vert = np.array([coord1, coord2, coord3])
        
        self._vertices.append(vert)
        
        
    def get_triangle(self, coord1, coord2, coord3):
        
        
        vert = np.array([coord1, coord2, coord3])
        
        #print self._triangles[0],'1'
        
        for triangle in self._triangles:
            

             
             if np.all(triangle._coord == vert):
        
        
                 return triangle 
    
    
        
    def create_triangle(self,  coord1, coord2, coord3):
           
            triangle = Triangle_container()
            #print triangle._node_id, 'new', triangle._data
            #print '1'
            #print triangle._data
            #triangle.add_data(data)
            
            vert = np.array([coord1, coord2, coord3])
            
            triangle.set_coord(vert)
            
            #triangle.add_node_id
            
            
            
            return triangle
            
            
            
        
    def add_triangle(self, coord1, coord2, coord3, triangle):
        
            vert = (coord1,coord2,coord3)
            
            self._triangles[vert] = triangle
            

        
        
        
            
            
         
#@profile
def pre_processing(i, num_data):
    
    start = time.time()
    

    grid = Triangle_table()
    
    n = 2**i+1
    
    int_n = n-2
    
    h =1/float(n-1)
    
    
    x = np.linspace(0, 1.0,num_data)
    y = np.linspace(0, 1.0,num_data)
    X, Y = np.meshgrid(x,y)
    
    model_func =Exy(X,Y)
    
    
    value = model_func.flatten()
    coordx = X.flatten()
    coordy = Y.flatten()
    Coord = zip(coordx, coordy)
    
    count =0
    for i in range(len(Coord)):

 
        data = Coord[i]
        #print data
        data_value = value[i]
        
        data_3D = data+(data_value,)
        
        data_3D = np.array([data_3D])
    
        if (data[0] !=0) and (data[0]!=1) and (data[1]!= 0 ) and (data[1] != 1): 
            #print data, 'data'
            
            # from 0 to n-2
            num_x = int(data[0]/float(h))
    
            num_y = int(data[1]/float(h))
        
        
            left_upper = (num_x*h, (num_y+1)*h)
            
            left_upper_id = num_x*n + (num_y+1)
            #print data, num_x, num_y, n,  left_upper_id
        
           
      
       
            left_down = (num_x*h, num_y*h)
            
            left_down_id = num_x*n + num_y
            
 
            
            
            #print left_down, ld_id
        
            right_upper = ((num_x+1)*h, (num_y+1)*h)
            
            right_upper_id = (num_x+1)*n+(num_y+1)
            
          
            #right_upper_int_id = num_x*int_n+num_y
            
   
        
            right_down = ((num_x+1)*h, num_y*h)
            
            right_down_id = (num_x+1)*n +  num_y
            
            
            
            
            if   (num_x != n-2)   and (num_y != n-2):
                
                right_upper_int_id = num_x*int_n+num_y
                
                
            else:
                
                right_upper_int_id = -1
                
               
                
                
            if (num_x != 0) and (num_y !=0):
                
                left_down_int_id = (num_x-1) * int_n +(num_y-1)
                
            else:
                
                left_down_int_id = -1
                
                
            if (num_x != n-2) and (num_y != 0):
                
                right_down_int_id = (num_x)*int_n + (num_y -1)
                
        
            else:
                
                right_down_int_id = -1
                
                
            if (num_x != 0) and (num_y != n-2):
                
                left_upper_int_id =  left_upper_int_id = (num_x-1) *int_n + num_y
                
            else:
                
                left_upper_int_id =-1
                
                
                
                
            
            # If the data is in the lower triangle
            if In_triangle(left_down, right_down, right_upper, data):
                
                #print left_down, right_down, right_upper
                
#                left_upper_id = num_x*n + (num_y+1)
#                
#                left_down_id = num_x*n + num_y
#                
#                right_upper_id = (num_x+1)*n+(num_y+1)
#                
#                right_down_id = (num_x+1)*n +  num_y
#            
      
            
           
            
                vert = (left_down, right_down, right_upper)
                #node_id = (left_down_id, right_down_id, right_upper_id)
                
                # If the triangle has been created and sotred
                try:
                    
                  
                    
                    tri = grid._triangles[vert]
                    
#                    print vert 
#                    
                    #print left_down_id, right_down_id,right_upper_id
                        
                    tri.add_data(data_3D)
                   
                    tri.add_node_id(left_down_id)
                    
                    tri.add_node_id(right_down_id)
                    
                    tri.add_node_id(right_upper_id)
                    
                    
                    tri.add_int_node_id(left_down_int_id)
                    
                    tri.add_int_node_id(right_down_int_id)
                    
                    tri.add_int_node_id(right_upper_int_id)
                    
                    
                    
                    
                    
                    
                # For a new triangle
                except KeyError:
                    
                    #print vert, '1'
                    
                    
                 
                    tri = grid.create_triangle(left_down, right_down, right_upper)
                    #print left_down_id, right_down_id,right_upper_id, right_down_int_id
                   
                    tri.add_data(data_3D)
                    
                    tri.add_node_id(left_down_id)
                    
                    tri.add_node_id(right_down_id)
                    
                    tri.add_node_id(right_upper_id)
                    
                    tri.add_int_node_id(left_down_int_id)
                    
                    tri.add_int_node_id(right_down_int_id)
                    
                    tri.add_int_node_id(right_upper_int_id)
                    
                    #print tri._node_id, 'id'
                
                    grid.add_triangle(left_down, right_down, right_upper,tri)
                    
                   
            

            # If the data is in the upper triangle 
            #else:
            if In_triangle(left_down, left_upper, right_upper, data):
                
                vert = (left_down, left_upper, right_upper)
                #print vert
                try:
                    
                    tri = grid._triangles[vert]
                    
                    #print left_down_id, left_upper_id,right_upper_id, left_upper_int_id
                    
                    tri.add_data(data_3D)
                    
                    tri.add_node_id(left_down_id)
                    
                    tri.add_node_id(left_upper_id)
                    
                    tri.add_node_id(right_upper_id)
                    
                    tri.add_int_node_id(left_down_int_id)
                    
                    tri.add_int_node_id(left_upper_int_id)
                    
                    tri.add_int_node_id(right_upper_int_id)
                

                
                    
                except KeyError:
                    #print vert, '1'

                 
                    tri = grid.create_triangle(left_down, left_upper, right_upper)
                    #print left_down_id, left_upper_id,right_upper_id, left_upper_int_id
                
                    tri.add_data(data_3D)
                    
                    tri.add_node_id(left_down_id)
                    
                    tri.add_node_id(left_upper_id)
                    
                    tri.add_node_id(right_upper_id)
                    
                    tri.add_int_node_id(left_down_int_id)
                    
                    tri.add_int_node_id(left_upper_int_id)
                    
                    tri.add_int_node_id(right_upper_int_id)
                    
                    grid.add_triangle(left_down,  left_upper, right_upper,tri)
    

    dictlist = []              
    for key ,value in grid._triangles.iteritems():
        
        temp = [key,value]
        
        dictlist.append(temp)

    done = time.time()
    elapsed = done - start
    #print(elapsed)
    return dictlist
            
        
        
        
        
    



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

    

def Lmatrix(i, num_data):
    
    grid = pre_processing(i, num_data)
    n = 2**i +1
    
    int_n = n-2
    size_L = int_n**2
    Lmatrix = lil_matrix((int_n**2, int_n**2))
    
    count = 0
    for tri in grid:

               count+=1
               #print count
               # node 1 should be the right angle node
               node1 = tri[0][0]
               idi = Coord_to_Id(node1,i)
               
               node2 = tri[0][1]
               idj = Coord_to_Id(node2,i)
              
               node3 = tri[0][2]       
               idk = Coord_to_Id(node3,i)
#               print (node1,node2, node3), (idi, idj, idk), 'tri'
               
               if 0<=idi<= size_L and 0<= idj<=size_L and 0<=idk<=size_L:
                   print (node1,node2, node3), (idi, idj, idk), 'tri'
                   Lmatrix[idi, idi] =4
                   Lmatrix[idj, idj] =4
                   Lmatrix[idk, idk] =4
               
                   if abs(idi-idj) ==1:
                       print idi, idj, 'ij'
                       
                       Lmatrix[idi, idj] =-1
                       Lmatrix[idj, idi] = -1 
                   

                       
                   if abs(idj-idk) ==1 :
                       print idj,idk, 'jk'
                       Lmatrix[idj,idk] =-1
                       Lmatrix[idk,idj] =-1
                       
                   if abs(idi-idj) == int_n:
                       
                       Lmatrix[idi, idj] =-1
                       Lmatrix[idj,idi] =-1
                       
                   if abs(idj-idk) ==int_n:
                        
                        Lmatrix[idj,idk] = -1
                        Lmatrix[idk, idj] = -1
               
    return Lmatrix.todense()
               
               

            
            
        
        
def G1matrix(i, num_data):
    
    grid = pre_processing(i, num_data)
    n = 2**i +1
    h =1/float(n-1)
    int_n = n-2
    size_L = int_n**2
    G1matrix = lil_matrix((int_n**2, int_n**2))
    
    count = 0
    for tri in grid:

               count+=1
               #print count
               # node 1 should be the right angle node
               node1 = tri[0][0]
               idi = Coord_to_Id(node1,i)
               
               node2 = tri[0][1]
               idj = Coord_to_Id(node2,i)
              
               node3 = tri[0][2]       
               idk = Coord_to_Id(node3,i)
#               print (node1,node2, node3), (idi, idj, idk), 'tri'
               
               if 0<=idi<= size_L and 0<= idj<=size_L and 0<=idk<=size_L:
                   print (node1,node2, node3), (idi, idj, idk), 'tri'
#                   G1matrix[idi, idi] =0*h
#                   G1matrix[idj, idj] =
#                   G1matrix[idk, idk] =6*h
               
                   if abs(idi-idj) ==1:
                       print idi, idj, 'ij'
                       
                       G1matrix[idi, idj] =-1*h
                       G1matrix[idj, idi] = 1*h
                   

                       
                   if abs(idj-idk) ==1 :
                       print idj,idk, 'jk'
                       G1matrix[idj,idk] =-1*h
                       G1matrix[idk,idj] =1*h
                       
                   if abs(idi-idj) == int_n:
                       
                       G1matrix[idi, idj] =2*h
                       G1matrix[idj,idi] =-2*h
                       
                   if abs(idj-idk) ==int_n:
                        
                        G1matrix[idj,idk] = 2*h
                        G1matrix[idk, idj] = -2*h
                        
                   if abs(idi - idk) == int_n+1:
                       
                       G1matrix[idi, idk] = 1*h
                       G1matrix[idk, idi] = -1*h
               
    return G1matrix.todense()/float(6)




    
    
    
#
#if __name__ == '__main__':
#    pre_processing(9,700)
#        
        
