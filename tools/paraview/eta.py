#!/usr/bin/env pvpython
r"Script for calculating viscosity coefficient from Couette flow problem between two parallel plates"

import sys
if len(sys.argv) != 2:
	raise NameError, 'Illegal number of arguments'
import os.path
if not os.path.isfile(sys.argv[1]):
	raise NameError, 'File ' + sys.argv[1] + ' doesn\'t exist' 

from paraview.simple import *
from numpy import *

all_cells = OpenDataFile(sys.argv[1])
data = servermanager.Fetch(all_cells)
L = data.GetBounds()[1]
N = data.GetNumberOfCells()
h = L/N
half_cells = ExtractCellsByRegion(all_cells)
half_cells.IntersectWith.Origin[0] = (L+h)/2
mean = servermanager.Fetch(IntegrateVariables(half_cells)).GetCellData()


flow = mean.GetArray('mass_flow').GetTuple3(0)[2]/h**3/(N/2)
qflow = mean.GetArray('heat_flow').GetTuple3(0)[2]/h**3/(N/2)
shear = data.GetCellData().GetArray('shear_stress').GetTuple3(0)[1]

U=.01

#print 'L = %.4e' % L
#print 'N = %.4e' % N
print 'flow = %.4e' % (flow/U)
print 'qflow = %.4e' % (qflow/U)
print 'shear = %.4e' % (shear/U)
