#!/usr/bin/env python
r"Plotting 2D distribution function cross-section"

### input parameters for script

axis = 2							# 0 -- X, 1 -- Y, 2 -- Z
value = 0							# cross section value

### loadind data from file

import sys
if len(sys.argv) != 2:						# program should have 1 argument = filename
        raise NameError, 'Illegal number of arguments'
import os.path
if not os.path.isfile(sys.argv[1]):				# check if file exist
        raise NameError, 'File ' + sys.argv[1] + ' doesn\'t exist'
import numpy as np
data = np.loadtxt(sys.argv[1], skiprows=4)			# get distribution function from file


### creating 2D projection

ax1, ax2 = (axis+1)%3, (axis+2)%3				# other axis for projection
f = open(sys.argv[1])
line = f.readline()
R, cut = line.split()					# load parameters from first line in the file
R, cut = int(R), float (cut)
line = f.readline()
dens, temp = map(float, line.split())
line = f.readline()
flow = map(float, line.split())
f.close()
delta = cut/R							# vel_grid step

def i2xi(i):
	return (i+.5-R)*delta

vel = np.fromfunction(lambda i: i2xi(i), (2*R,))		# 1D vel_grid
value = vel[np.abs(vel-value).argmin()]				# the nearest vel_grid node to the specific value
f = np.zeros((2*R, 2*R))					# 2D projection of vel_grid

for l in data[1:]:
	xi = round(l[ax1]/delta-.5+R)				# first coordinate
	yi = round(l[ax2]/delta-.5+R)				# second coordinate
	if np.abs(l[axis]-value)*1e4 < delta:
		f[xi,yi] = l[3]					# filling array

f = np.ma.masked_where(f == 0, f)				# masking nodes not included in the vel_grid ball

maxw = dens/(np.pi*temp)**1.5 * np.fromfunction(lambda i,j: \
	np.exp(-((i2xi(i)-flow[ax1])**2+(i2xi(j)-flow[ax2])**2+(value-flow[axis])**2)/temp), (2*R, 2*R))

### plotting 3D graph

print dens, temp
print flow[ax1]

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x, y = np.outer(vel,np.ones(2*R)), np.outer(np.ones(2*R),vel)	# filling coordinates for function "plot"
ax.plot_wireframe(x, y, f-maxw, rstride=1, cstride=1, color='b')
#ax.plot_surface(x, y, f, rstride=1, cstride=1, color='w')

ax.set_xlabel(chr(ord('X') + ax1))
ax.set_ylabel(chr(ord('X') + ax2))
ax.set_zlabel(chr(ord('X') + axis)+' slice= '+str(value))
plt.show()

