#!/usr/bin/env pvpython
r"Script for calculating thermal conductivity coefficient from heat transfer problem between two parallel plates"

import sys
if len(sys.argv) != 2:
	raise NameError, 'Illegal number of arguments'
import os.path
if not os.path.isfile(sys.argv[1]):
	raise NameError, 'File ' + sys.argv[1] + ' doesn\'t exist' 

from paraview.simple import *
from numpy import *

all_cells = OpenDataFile(sys.argv[1])
outline = servermanager.Fetch(Outline(all_cells))
L = outline.GetBounds()[1]

radius = L/10
center = L/2
line = PlotOverLine(all_cells)
line.Source.Point1[0] = center - radius
line.Source.Point2[0] = center + radius
points = 10
line.Source.Resolution = 2*points
data = servermanager.Fetch(line).GetPointData()

dens = data.GetArray('density').GetTuple1(points)
temp = data.GetArray('temperature').GetTuple1(points)
flow = data.GetArray('mass_flow').GetTuple3(points)[0]
qflow = data.GetArray('heat_flow').GetTuple3(points)[0]
grad_temp = (data.GetArray('temperature').GetTuple1(2*points) - data.GetArray('temperature').GetTuple1(0))/(2*radius)

dT=.01
kappa = -qflow/grad_temp/sqrt(temp)
kappa_all = -qflow/dT/sqrt(temp) * L * (1+2.4001*sqrt(pi)/L)
heat_flow = -(2+dT)/dT*qflow/(dens*temp**1.5)

print 'L = %d' % L
print 'density = %.4e' % dens
print 'Sone_qflow = %.4e' % heat_flow
#print 'temperature = %.2e' % temp
print 'grad_temp = %.2e' % grad_temp
print 'heat_flow = %.2e' % qflow
print 'kappa_all = %.4f' % round(kappa_all, 4)
print 'mass_flow = %.1e' % flow
print 'kappa = %.4f' % round(kappa, 4)
