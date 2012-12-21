#!/usr/bin/env pvpython
r"Script for calculating plane Poiseuille mass flow in the center cross-section"

import sys
if len(sys.argv) != 2:
	raise NameError, 'Illegal number of arguments'
import os.path
if not os.path.isfile(sys.argv[1]):
	raise NameError, 'File ' + sys.argv[1] + ' doesn\'t exist' 

from paraview.simple import OpenDataFile, PlotOverLine, Outline, servermanager
all_cells = OpenDataFile(sys.argv[1])
outline = servermanager.Fetch(Outline(all_cells))
L = outline.GetBounds()[1]
d = 2*outline.GetBounds()[5]

radius = 5
center = L/2 + .5
line = PlotOverLine(all_cells)
line.Source.Point1[0] = center - radius
line.Source.Point2[0] = center + radius
points = 20*radius
line.Source.Resolution = points

sum_flow = 0
aver_dens = 0
for i in range(d/2):
	line.Source.Point1[2] = line.Source.Point2[2] = i + .5
	data = servermanager.Fetch(line).GetPointData()
	flow = data.GetArray('mass_flow').GetTuple3(points/2)[0]
	dens = data.GetArray('density').GetTuple1(points/2)
	press = sum(data.GetArray('pressure').GetTuple3(points/2))/3
	grad_press = (sum(data.GetArray('pressure').GetTuple3(points)) 
	- sum(data.GetArray('pressure').GetTuple3(0)))/3/2/radius
	sum_flow += flow * press / grad_press / dens
	aver_dens += dens

aver_dens = aver_dens/(d/2)
Q = -2**.5 * sum_flow / d**2

print 'd = %d' % d
print 'L = %d' % L
print 'rho_av = %.4f' % round(aver_dens, 4)
print 'Q = %.4f' % round(Q, 4)
