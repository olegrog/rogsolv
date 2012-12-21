#!/usr/bin/env python
r"Plotting 2D distribution function cross-section"

import sys
if len(sys.argv) != 2:
        raise NameError, 'Illegal number of arguments'
import os.path
if not os.path.isfile(sys.argv[1]):
        raise NameError, 'File ' + sys.argv[1] + ' doesn\'t exist'

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

R, cut = 16, 3.4
delta = cut/R
dim = (2*R, 2*R)
vel = np.fromfunction(lambda i: (i+.5-R)*delta, (2*R,))
x, y = np.outer(np.ones(2*R),vel), np.outer(vel,np.ones(2*R))
f = np.zeros(dim)
f[:,:] = 0

data = np.genfromtxt(sys.argv[1])
for l in data:
	xi = round(l[0]/delta-.5+R)
	yi = round(l[1]/delta-.5+R)
	f[xi,yi] = l[2]

#X = data[:,0]
#Y = data[:,1]
#Z = data[:,2]
#X, Y, Z = axes3d.get_test_data(0.1)

f = np.ma.masked_where(f == 0, f)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, f, rstride=1, cstride=1, color='b')
#ax.plot_surface(x, y, f, rstride=1, cstride=1, color='w')

plt.show()

