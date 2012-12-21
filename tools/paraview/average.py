#!/usr/bin/env pvpython
r"Script for calculating time average macroparameters"

import sys
if len(sys.argv) != 3:
    raise NameError, 'Illegal number of arguments. Parameters: path field. Example: average.sh ./minmod/ mass_flow'
import os.path
directory = sys.argv[1]
if not os.path.isdir(directory):
    raise NameError, 'File ' + directory + ' doesn\'t exist' 

from paraview.simple import *
import numpy

schemes = ["first" "mc"]
flag = True
time_range = range(1800,2001,10)
for time in time_range:
    all_cells = OpenDataFile(directory+"/VTK/"+str(time)+"data.vtk")
    data = servermanager.Fetch(all_cells)
    array = data.GetCellData().GetArray(sys.argv[2])
    if flag:
        L = array.GetNumberOfTuples()
        a = numpy.zeros(L, dtype = numpy.float128)

    flag = False
    for i in range(L):
        a[i] += array.GetTuple3(i)[0]

for i in range(L):
    a[i] /= len(time_range)
    print a[i]
