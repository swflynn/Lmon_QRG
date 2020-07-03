#=============================================================================80
#                     2D Contour Plot for Water Potentials
#=============================================================================80
#   Discussion:
#Given a set of x1,x2,V data compute a contour plot
#==============================================================================#
#   Modified:
#20 June 2020
#   Author:
#Shane Flynn
#==============================================================================#
#Make it compute the contour (pass in data file from uniform grid)
#then pass in data file for any grid and plot the grid itself over the contour
#==============================================================================#
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
#==============================================================================#
N=10000

x, y, z = np.genfromtxt(r'cartesian_rU.dat', unpack=True)
xll = x.min();  xul = x.max();  yll = y.min();  yul = y.max()

xi = np.linspace(xll, xul, N)
yi = np.linspace(yll, yul, N)
zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

contours = plt.contour(xi, yi, zi, 6, colors='black')
#plt.show()
plt.savefig('cartesian_rU.png')
