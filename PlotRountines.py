#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:50:15 2019

@author: shilu
"""

def smoother_plot(u):
    """ Plot the approximation of u in 3D. This routine can be used in
    pre-smooth and post-smooth stages"""
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    #from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np
    from functions import Linear, Xlinear, Ylinear, Zero
    
    # Get number of nodes
    [xdim, ydim] = u.shape
    
    fig = plt.figure(figsize=(8,5))
    ax = fig.gca(projection='3d')
    
    
    # Get the spacing, for example, 5 nodes on side yields the spacing h =1/4
    h=1/float(xdim-1)
    
    x1, y1 = np.meshgrid(np.arange(h, 1, h), np.arange(h, 1, h))
    
#    true_u = Xlinear(x1,y1)
    
#    error = u[0][1:-1,1:-1]- true_u
    ax.plot_surface(x1, y1, u[1:-1,1:-1] ,cmap='viridis',linewidth=0)

    
    # Make the ticks looks pretty
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    plt.show()    
    
    
    
def convergence_plot(cyclenumber,norm):
    """ Plot the convergence rate for V-cycle. We have options for error convergence
    and residue convergence
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from pylab import title, xlabel, ylabel, clf, plot,show, legend
    
    
    xline = np.arange(cyclenumber+1)
    plt.figure(figsize=(4,5))
    plt.semilogy(xline, norm, 'bo-', xline, norm, 'k',label='sdad')
    #plt.semilogy(xline, egg1norm, 'bo', xline, egg1norm, 'k',label='sdad')
    title('Convergence with Residual(Richardson)')
    xlabel('Number of cycles')
    ylabel('Error under l2 norm')
    plt.show()
    
    
    
#fig,ax = subplots()
#
#ax.scatter(Eigenvalue5_D.real,Eigenvalue5_D.imag, s=1)
#ax.set_xlim([-1, 1])
#ax.set_ylim([-10, 10])