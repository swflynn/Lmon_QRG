#=============================================================================80
#                     2D Contour Plot for Water Potentials
#=============================================================================80
#       Discussion:
#Pass in 2 data files, data for potential, data for grid
#==============================================================================#
#To do:
#Read in data to make Potential
#Plot Ecut contour on this poitential
#Plot a data file over the potential
#==============================================================================#
#       Modified:
# 26 June 2020
#       Author:
# Shane Flynn
#==============================================================================#
import sys
import pandas as pd
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
#==============================================================================#
#script             ==>Name of this script
#pot_data           ==>Data file with x y V evaluations
#grid_data          ==>Data file with x y evaluations to plot over potential
#x/y min/max        ==>Range from pot_data
#num                ==>Number of points to interpolate (xi*yi)
#xi/yi              ==>points at which to interpolate
#==============================================================================#
#                       Read in Potential/Grid Files
#==============================================================================#
script,pot_data,grid_data=sys.argv
num=int(input())
levels=int(input())
df_pot=pd.read_csv(pot_data,delim_whitespace=True,header=None,dtype=np.float64)
df_grid=pd.read_csv(grid_data,delim_whitespace=True,header=None,dtype=np.float64)
x=df_pot[0]; xmin=df_pot[0].min(); xmax=df_pot[0].max()
y=df_pot[1]; ymin=df_pot[1].min(); ymax=df_pot[1].max()
z=df_pot[2]
#==============================================================================#
#                          Interpolation Points
#==============================================================================#
xi=np.linspace(xmin,xmax,num)
yi=np.linspace(ymin,ymax,num)
zi=scipy.interpolate.griddata((x, y),z,(xi[None,:],yi[:,None]),method='cubic')
#==============================================================================#
#                                grid_data
#==============================================================================#
x_grid=df_grid[0]
y_grid=df_grid[1]
#==============================================================================#
#                             Plot Grid + Contours
#==============================================================================#
autokcalmol=627.5096
#level2=[1/autokcalmol,5/autokcalmol,10/autokcalmol,20/autokcalmol]
level2=[10/autokcalmol]
#==============================================================================#
contours=plt.contour(xi,yi,zi,levels,colors='black')
plt.contour(xi,yi,zi,level2,colors='red') #quick fix to plot Ecut
plt.plot(x_grid,y_grid,'bo')
plt.show()
###plt.savefig(png_file)
#==============================================================================#
#                                Output
#==============================================================================#
f=open('output','w')
f.write('2D water potential contour plot\n')
f.write(f'Input potential file: {pot_data}\n')
f.write(f'Input grid file: {grid_data}\n')
f.write(f'number of interpolation evaluations (1d): {num}')
f.write(f'xmin: {xmin}, xmax: {xmax}, ymin: {ymin}, ymax: {ymax}')
f.close()
