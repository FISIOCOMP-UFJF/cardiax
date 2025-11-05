#!/usr/bin/env python

import sys, gzip, pdb, time
import numpy as np
from sctools import read_array_pts, read_array_elem
from sctools import carp_map_label_node as carp_map

def carp2vtk (carpMesh, outputMesh):

    ptsFile  = carpMesh + '.pts'
    elemFile = carpMesh + '.elem'
    
    pts  = read_array_pts (ptsFile)
    elem = read_array_elem (elemFile)

    elem = np.array(elem)
    elconec = np.array(elem[:,1:-1], dtype=np.int)
    eltypes = elem[:,0]
    elregio = np.array(elem[:,-1], dtype=np.int)
   
    numPts  = np.shape(pts)[0]
    numElem = np.shape(elem)[0]

    out_vtk = file(outputMesh, 'w')

    #out_vtk.write('# vtk DataFile Version 2.0\n')
    
    out_vtk.write('# vtk DataFile Version 1.0\n')
    out_vtk.write('Auto generated VTK file from input carp mesh (carp2vtk)\n')
    out_vtk.write('ASCII\n\n')

    out_vtk.write('DATASET UNSTRUCTURED_GRID\n')

    out_vtk.write('POINTS %d float\n' % numPts)
    for i in xrange(numPts):
        out_vtk.write('%f %f %f\n' % (pts[i,0],pts[i,1],pts[i,2]))
    out_vtk.write('\n')

    cell_size = get_cell_size(numElem, eltypes)
    
    out_vtk.write('CELLS %d %d\n' % (numElem,cell_size))
    for i in xrange(numElem):
        typestr = str(eltypes[i]).lstrip()
        elmtype = carp_map[typestr]
        out_vtk.write('%d ' % elmtype)
        for node in xrange(elmtype):
            out_vtk.write('%d ' % elconec[i,node])
        out_vtk.write('\n')
    out_vtk.write('\n')

    out_vtk.write('CELL_TYPES %d\n' % numElem)
    for i in xrange(numElem):
        out_vtk.write('%d\n' % 5)
#        out_vtk.write('%d\n' % elregio[i])
    out_vtk.write('\n')

    out_vtk.write('POINT_DATA %d\n' % numPts)
    out_vtk.write('SCALARS scalars float\n')
    out_vtk.write('LOOKUP_TABLE default\n')
    for i in xrange(numPts):
        out_vtk.write('%f\n' % 1.0)
    out_vtk.write('\n')

    out_vtk.close()

def get_cell_size(num_elem, el_types):    
    sum = 0
    for i in xrange(num_elem):
        typestr = str(el_types[i]).lstrip()
        sum += carp_map[typestr]

    sum = sum + num_elem
    return sum   

if __name__ == "__main__":

    np.set_printoptions(precision=5)

    if len(sys.argv) < 3:
        print( "\n Usage: carp2vtk <carp_mesh> <vtk_mesh>\n")
        sys.exit(1)

    carp2vtk(sys.argv[1], sys.argv[2])


# end of main
