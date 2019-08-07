#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 18:07:40 2018

@author: shilu

The test is aiming to approximate a two dimensional surface


"""
import numpy as np
from numpy import  pi, sin, cos, exp, inf
from scipy import zeros, linalg



from grid.Grid import Grid
from BuildSquare import build_square_grid, build_square_grid_matrix
from BuildEquation import build_equation_linear_2D, set_polynomial_linear_2D,\
Poisson_tri_integrate, TPS_tri_intergrateX, TPS_tri_intergrateY, NotIn_triangle, NIn_triangle, build_matrix_fem_2D

from Triangle import triangle_iterator, Interior_triangle_iterator
from grid.NodeTable import node_iterator, not_slave_node
from grid.ConnectTable import connect_iterator 

from grid.function.FunctionStore import zero,exp_soln, linear, plain
import matplotlib.pyplot as plt
from matplotlib.mlab import bivariate_normal
from mpl_toolkits.mplot3d import Axes3D

from PlotPackman import plot_fem_grid
from scipy.sparse.linalg import spsolve

from copy import copy
from operator import itemgetter


def sin_soln(x,y):
    """Define the function sin(pi x_0)sin(pi x_1)"""
    return sin(pi*x)*sin(pi*y)

def Plain(x,y):
    
    return x-x+1

def l2(x,y):
    
    return x**2 + y**2

def Linear(x,y):
    
    return x+y


def expn(x,y):
    """Define the function exp(x_0)exp(x_1)"""
    return exp(x)*exp(y)
    



#Create grid and multivariate normal
x = np.linspace(0, 1.0,30)
y = np.linspace(0, 1.0,30)

X, Y = np.meshgrid(x,y)

Z2 = Linear(X,Y)

#Make a 3D plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z2,cmap='viridis',linewidth=0)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()
#
data = Z2.flatten()
coordx = X.flatten()
coordy = Y.flatten()
Coord = zip(coordx, coordy)
Coord1 = copy(Coord)



##############################################################################
##############################################################################

grid = Grid()

# Build a 5*5 grid
i = 2
n = 2**i

alpha =1

true_soln = plain

build_square_grid(n, grid, true_soln)
#build_square_grid(n, grid, true_soln)





Nl=[]
for node in (not_slave_node(grid)):
    #print(i)
    Nl.append([node,node.get_node_id()])
Nl=sorted(Nl, key = itemgetter(1))
Nl=[node[0] for node in Nl]
        

        







def generate_Amatrix():
    
    
    Nl=[]
    for node in (not_slave_node(grid)):
        
        Nl.append([node,node.get_node_id()])
    Nl=sorted(Nl, key = itemgetter(1))
    Nl=[node[0] for node in Nl]
            
    
    for node in Nl:
        node.set_value(Nl.index(node))
    
    Amatrix = zeros([len(Nl), len(Nl)])
    
    #plot_fem_grid(grid)
    
    for node in node_iterator(grid):
         if not node.get_slave():
             print node.get_node_id()._id_no, node.get_value()

    for tri in Interior_triangle_iterator(grid):
        
    
        
        if (tri[1].get_node_id() < tri[2].get_node_id()):
            
            #print tri[0].get_node_id()._id_no, tri[1].get_node_id()._id_no, tri[2].get_node_id()._id_no
            
                            
            basis1 = set_polynomial_linear_2D(tri[0], tri[2], tri[1])
                    
            basis2 = set_polynomial_linear_2D(tri[1], tri[2], tri[0])
                    
            basis3 = set_polynomial_linear_2D(tri[2], tri[1], tri[0])
            

            
            for i in Coord:
                

            
                if  NIn_triangle(tri[0] , tri[1], tri[2], i):
                    
                    #print In_triangle(tri[0] , tri[1], tri[2], i)
                    
                    #print tri[0].get_node_id()._id_no, tri[1].get_node_id()._id_no, tri[2].get_node_id()._id_no, i,\
                    #tri[0].get_coord(), tri[1].get_coord(), tri[2].get_coord()
                
                    
                    
                    
                    Idi = int(tri[0].get_value())
                    
                    
                    Amatrix[Idi, Idi] += basis1.eval(i[0], i[1]) * basis1.eval(i[0], i[1])
            

                    if not tri[1].get_slave():
                        
                        #print tri[1].get_node_id()._id_no, tri[1].get_value()
                        
                        #print tri[1].get_value()

                        Idj = int(tri[1].get_value())
                        
                        

                        
                        Amatrix[Idi,Idj] += basis1.eval(i[0], i[1]) * basis2.eval(i[0], i[1])
                        
                
                    if not tri[2].get_slave():
                        
                        #print tri[2].get_node_id()._id_no
                        
                        Idk = int(tri[2].get_value())
                        
                        
                        
                        Amatrix[Idi, Idk] += basis1.eval(i[0], i[1]) * basis3.eval(i[0], i[1])
                        
                    # Remove the data once it's been used 
                    #Coord.remove(i)
            
#        
    
    return Amatrix/float(len(Coord))

Amatrix = generate_Amatrix()

    
    
    
    
##############################################################################
##############################################################################    

def generate_dmatrix():
    
    """ This rountine gernerates dmatrix after deducted boundary conditions
    """


    Nl=[]
    for node in (not_slave_node(grid)):
    
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
            
    
    for node in Nl:
        node.set_value(Nl.index(node))
    
    dmatrix = zeros([len(Nl),1])
    
    #print Amatrix

    for tri in Interior_triangle_iterator(grid):
        
        if (tri[1].get_node_id() < tri[2].get_node_id()):
            
            #print tri[0]
            
            #print tri[0].get_node_id()._id_no, tri[1].get_node_id()._id_no, tri[2].get_node_id()._id_no
            basis1 = set_polynomial_linear_2D(tri[0], tri[2], tri[1])
            
            Idi = int(tri[0].get_value())
            
            for i in range(len(Coord)):
                
                
                if  NIn_triangle(tri[0] , tri[1], tri[2], Coord[i]):
                    
                    #print  basis1.eval(Coord[i][0], Coord[i][1]), data[i]
                
                    dmatrix[Idi,0] += basis1.eval(Coord[i][0], Coord[i][1]) * data[i]
                    #print dmatrix[Idi,0]
                    
    #print dmatrix/float(len(Coord))
                    
                        
    return dmatrix/float(len(Coord))
            

dmatrix1 =  generate_dmatrix()
dmatrix =    dmatrix1-np.array([[0.0275862],
                               [0.0344828],
                               [0.0344828],
                               [0.0275862]])





##############################################################################
##############################################################################






def generate_Lmatrix():
    build_equation_linear_2D(grid, Poisson_tri_integrate, zero)
    #build_equation_linear_2D(grid, Poisson_tri_integrate, zero)

    Nl=[]
    for node in (not_slave_node(grid)):
        #print(i)
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
            
    
    for node in Nl:
        node.set_value(Nl.index(node))
        
        
    Lmatrix  = zeros([len(Nl),len(Nl)])

    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    stiffness = grid.get_matrix_value(node.get_node_id(), endpt1)[0] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        Lmatrix[i, j] = stiffness
    return Lmatrix
            
  
Lmatrix = generate_Lmatrix() 
#Lmatrix = np.array([[ 4., -1., -1.,  0.],
# [-1.,  4.,  0., -1.],
# [-1.,  0.,  4., -1.],
# [ 0., -1., -1.,  4.]])
        
        


##############################################################################
##############################################################################
def generate_XGmatrix():
    
    build_equation_linear_2D(grid, TPS_tri_intergrateX, zero)


        
    
    Nl=[]
    for node in (not_slave_node(grid)):
        #print(i)
        Nl.append([node,node.get_node_id()])
        
    Nl=sorted(Nl, key = itemgetter(1))
    
    Nl=[node[0] for node in Nl]
        
    
    for node in Nl:
        node.set_value(Nl.index(node))
        
    XGmatrix  = zeros([len(Nl),len(Nl)])

    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    stiffness = grid.get_matrix_value(node.get_node_id(), endpt1)[0]
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        XGmatrix[i, j] = stiffness
    return XGmatrix


XGmatrix=generate_XGmatrix()
       
        


##############################################################################
##############################################################################        
        
def generate_YGmatrix():
    
    build_equation_linear_2D(grid, TPS_tri_intergrateY, zero)

    Nl=[]
    for node in (not_slave_node(grid)):
        #print(i)
        Nl.append([node,node.get_node_id()])
    Nl=sorted(Nl, key = itemgetter(1))
    Nl=[node[0] for node in Nl]
        
    
    for node in Nl:
        node.set_value(Nl.index(node))
        
    YGmatrix  = zeros([len(Nl),len(Nl)])

    for node in node_iterator(grid):
    
        # Ignore slave (or boundary) nodes
        if not node.get_slave():
            
            # Which row corresponds to the current node?
            i = int(node.get_value())
        
            for endpt1 in connect_iterator(grid, node.get_node_id()):
    
                    # Which column corresponds to the current node?
                    j = int(grid.get_value(endpt1))
                    
                    # What is the corresponding matrix value (in the FEM grid)
                    stiffness = grid.get_matrix_value(node.get_node_id(), endpt1)[0] 
        
                    # We must not include slave nodes in the matrix columns
                    if not grid.get_slave(endpt1):
                        YGmatrix[i, j] = stiffness
    return YGmatrix        
        
        
        
        
YGmatrix=generate_YGmatrix()    
    
##############################################################################
##############################################################################

ZeroMatrix = zeros([len(Nl),len(Nl)])



one = zeros([16,1])
one[0:4]=np.ones((4,1))

def Lets_Make_the_Damn_Big_Matrix():
    
    BigMat = np.block([[Amatrix, ZeroMatrix, ZeroMatrix, Lmatrix],\
                       [ZeroMatrix, alpha*Lmatrix, ZeroMatrix, XGmatrix.T],\
                       [ZeroMatrix, ZeroMatrix, alpha*Lmatrix,  YGmatrix.T],\
                       [Lmatrix, XGmatrix, YGmatrix, ZeroMatrix]])
    return BigMat
    
BigMat = Lets_Make_the_Damn_Big_Matrix()

dpline=zeros([4*len(Nl),1])
dpline[0:len(Nl)]=dmatrix

print dpline


#value = np.dot(linalg.inv(BigMat), dvector)
#print value
value_vector = linalg.solve(BigMat, dpline)[0:len(Nl)]
value_vector2 = linalg.solve(BigMat, dpline)
        




for node in node_iterator(grid):
    
    # The error is only calculated at the interior nodes
    
    if not node.get_slave():
        
        # What row in the matrix does the current node correspond to
        i = int(node.get_value())
        

        value = value_vector[i]
        #print value
        
        #print value

        
        # Record the value at the current node
        node.set_value(value)
        
        
        
    # If the node is a slave node, its value should be given by the
    # boundary condition
    else:
        node.set_value(plain(node.get_coord()))
        

    
  

def plot_fem_solution(grid):
    from grid.Grid import Grid
    from grid.NodeTable import node_iterator
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from matplotlib.pyplot import show, figure
    from scipy.interpolate import griddata
    from numpy import mgrid, array
#   

        
# The connection values are set separately with respect to the location, 

    #grid=Grid()
#    
    #Find the position of the nodes and the values
    node_x=[]
    node_y=[]
    node_v=[]
    for node in node_iterator(grid):
        coord = node.get_coord()
        node_x.append(coord[0])
        node_y.append(coord[1])
        node_v.append(node.get_value())
        
    # Store the results in an array
    
    node_x = array(node_x)
    node_y = array(node_y)
    node_v = array(node_v)
    
    #print('node_x',node_x)
    #print('node_value', node_v)
    
    
    # Initialise the figure
    fig = figure()
    ax = fig.gca(projection='3d') 
    ax = fig.gca()
    
    
    # Interpolate the nodes onto a structured mesh
    X, Y = mgrid[node_x.min():node_x.max():10j,
                 node_y.min():node_y.max():10j]
    
    Z = griddata((node_x,node_y), node_v, (X,Y), method='cubic')
    
    
    # Make a surface plot
    ax.plot_surface(X, Y, Z,cmap='viridis',linewidth=0)
    
    # Set the z axis limits
    ax.set_zlim(node_v.min(),node_v.max())
    
    # Make the ticks looks pretty
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    
    # Include a colour bar
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    
    # Show the plot
    show()

        
        
plot_fem_solution(grid)    
#    
