#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 12:21:39 2018

@author: shilu
"""


from numpy import  pi, sin, cos, exp

def sin_soln(x,y):
    """Define the function sin(pi x_0)sin(pi x_1)"""
    return sin(15*pi*x)*sin(15*pi*y)+sin(pi*x)*sin(pi*y)

def cos_sin(x, y):
    """Define the function pi cos(pi x_1)sin(pi y_1)"""
    
    return pi*cos(pi*x)*sin(pi*y)

def sin_cos(x, y):
    """Define the function pi sin(pi x_1)cos(pi y_1)"""
    
    return pi*sin(pi*x)*cos(pi*y)

def exp_soln(x,y):
    """Define the function exp(x_0)exp(x_1)"""
    return exp(x)*exp(y)
    

def exp_w(x,y):
    """Define the function -2*exp(x_0)exp(x_1)"""
    return -1.0*exp(x)*exp(y)


def Cos2(x,y):
    
    return cos(2*pi*x)*cos(2*pi*y)

def Xcos2(x,y):
    
    return -2*pi*cos(2*pi*y)*sin(2*pi*x)

def Ycos2(x,y):
    
    return -2*pi*cos(2*pi*x)*sin(2*pi*y)

def XYcos2(x,y):
    
    return (8*0.000000001*pi**2)*cos(2*pi*x)*cos(2*pi*y)



def Zero(x,y):
    
    return x-x

def Linear(x,y):
    
    return x+y

def Xlinear(x,y):
    
    return x-x+1

def Ylinear(x,y):
    
    return x-x+1
    


def Linear2x(x,y):
    
    return 2*x+y

def Xlinear2x(x,y):
    
    return x-x+2.0
    

def Ylinear2x(x,y):
    
    return y-y+1.0


def Cube(x,y):
    
    return x**3 + y**3

def Xcube(x,y):
    
    return 3*x**2


def Ycube(x,y):
    
    return 3*y**2 

def XYcube(x,y):
    
     return 6*(x+y)
 
    

def Exy(x,y):
    return exp(3/((x-0.5)**2+(y-0.5)**2+1))


def Xexy(x,y):
    return -6 * exp(3/((x-0.5)**2+(y-0.5)**2+1))*(-0.5+y)/((x-0.5)**2+(y-0.5)**2+1)**2
    
def Yexy(x,y):
    return -6 * exp(3/((x-0.5)**2+(y-0.5)**2+1))*(-0.5+x)/((x-0.5)**2+(y-0.5)**2+1)**2
    
def XYexy(x,y):
    
    return  1*(36 * exp(3/((x-0.5)**2+(y-0.5)**2+1))*(-0.5+x)**2/((x-0.5)**2+(y-0.5)**2+1)**4+ \
           24 * exp(3/((x-0.5)**2+(y-0.5)**2+1))*(-0.5+x)**2/((x-0.5)**2+(y-0.5)**2+1)**3-\
           12 * exp(3/((x-0.5)**2+(y-0.5)**2+1))/((x-0.5)**2+(y-0.5)**2+1)**2+\
           36 * exp(3/((x-0.5)**2+(y-0.5)**2+1))*(-0.5+y)**2/((x-0.5)**2+(y-0.5)**2+1)**4+\
           24 * exp(3/((x-0.5)**2+(y-0.5)**2+1))*(-0.5+y)**2/((x-0.5)**2+(y-0.5)**2+1)**3)
           
    