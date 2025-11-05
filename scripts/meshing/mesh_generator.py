import numpy as np
#import xml.etree.ElementTree as ET
import os, sys
import benchmark_ellipse_linear as bel



if (len(sys.argv) < 4):
    print("\n Usage: mesh_generator <filename> <nab> <ncirc> <ntrans>\n")
    sys.exit(-1)

filename = sys.argv[1]

h=[int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])]
print h
ninc=10
pressure=2
C = 1.1
bf = 6.6
bt = 4
bfs = 2.6
K=300
Ta = 0
top = 5.
outer_long_axis = 20.
outer_short_axis = 10.
wall_thickness = 3.
#fiber angles
fiber_angle_epi  = -60.
fiber_angle_endo = 60.

mesh_params = [top, outer_long_axis, outer_short_axis, wall_thickness, fiber_angle_epi, fiber_angle_endo]
model_params = np.array([C, bf, bt, bfs, Ta, K])

bel.createXml(filename, mesh_params, model_params, h, ninc, pressure)



