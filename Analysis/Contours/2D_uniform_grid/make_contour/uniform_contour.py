#=============================================================================80
#                     2D Contour Plot for Water Potentials
#=============================================================================80
#   Discussion:
#Given a set of x1,x2,V data compute a contour plot
#==============================================================================#
#   Modified:
#26 June 2020
#   Author:
#Shane Flynn
#==============================================================================#
#script        ==>Name of THIS script
#data_file     ==>Name of data file with coordiantes and energy
#png_file      ==>Name of the output png file
#N             ==>Number of points to construct the Contour
#==============================================================================#
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
#==============================================================================#
script,data_file,png_file=sys.argv
#==============================================================================#
N=input()
try:
    N=int(N)
except:
    print('N needs to be an integer, check input file')
    sys.exit(-1)
#==============================================================================#
x,y,z=np.genfromtxt(data_file, unpack=True)
xll=x.min(); xul=x.max(); yll=y.min(); yul=y.max()
#==============================================================================#
xi=np.linspace(xll,xul,N)
yi=np.linspace(yll,yul,N)
zi=scipy.interpolate.griddata((x, y),z,(xi[None,:],yi[:,None]),method='cubic')
#==============================================================================#
contours = plt.contour(xi, yi, zi, 6, colors='black')
plt.savefig(png_file)
#==============================================================================#
f=open('output','w')
f.write('2D water potential contour plot\n')
f.write(f'Input data file: {data_file}\n')
f.write(f'Output png file: {png_file}\n')
f.write(f'N: {N}')
f.close()
