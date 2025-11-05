#!/usr/bin/env python

import sys, time, os
import numpy as np
from sctools import read_array_pts, read_array_elem

def carp2Gmsh (carpMesh, outputMesh):

    ptsFile  = carpMesh + '.pts'
    elemFile = carpMesh + '.elem'
    
    pts  = read_array_pts (ptsFile)
    elem = read_array_elem (elemFile)

    numPts  = np.shape(pts)[0]
    numElem = np.shape(elem)[0]

    # basic header
    outputFile = file(outputMesh, 'w')
    outputFile.write('$MeshFormat\n2.0 0 8\n$EndMeshFormat\n')

    # nodes section
    outputFile.write('$Nodes\n')
    outputFile.write('%d\n' % numPts)
    for i in xrange(numPts):
        outputFile.write('%d %f %f %f\n' % (i+1,pts[i,0],pts[i,1],pts[i,2]))
    outputFile.write('$EndNodes\n')

    # elements section
    outputFile.write('$Elements\n')
    outputFile.write('%d\n' % numElem)
 
    for i in xrange(numElem):
        elem_type = str(elem[i][0])
        
        if (elem_type == "Tr"):
            write_tr(outputFile, i+1, np.array(elem[i][1:],dtype=int)+1 )
        elif (elem_type == "Tt"):
            write_tt(outputFile, i+1, np.array(elem[i][1:],dtype=int)+1 )
        elif (elem_type == "Qd"):
            write_qd(outputFile, i+1, np.array(elem[i][1:],dtype=int)+1 )
        elif (elem_type == "Hx"):
            write_hx(outputFile, i+1, np.array(elem[i][1:],dtype=int)+1 )
        else:
            print(" Error: elem_type not supported.")
            sys.exit(1)          

    outputFile.write('$EndElements\n')
    outputFile.close()

def write_qd(out, e, conec):
    out.write('%d 3 2 99 2 %d %d %d %d\n' % (e, 
        conec[0], conec[1], conec[2], conec[3]))

def write_tr(out, e, conec):
    out.write('%d 2 2 99 2 %d %d %d\n' % (e, 
        conec[0], conec[1], conec[2]))

def write_tt(out, e, conec):
    out.write('%d 4 2 99 2 %d %d %d %d\n' % (e, 
        conec[0], conec[1], conec[2], conec[3]))

def write_hx(out, e, conec):
    out.write('%d 5 2 99 2 %d %d %d %d %d %d %d %d\n' % (e, 
        conec[0], conec[1], conec[2], conec[3],
        conec[4], conec[5], conec[6], conec[7]))

if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: carp2gmsh <carp_mesh> <gmsh_output>\n")
        sys.exit(-1)

    # parse and check input
    carp_mesh = sys.argv[1]

    if (not os.path.isfile(carp_mesh+'.pts')) or (not os.path.isfile(carp_mesh+'.elem')):
        print("\n Error: the input carpfile %s does not exist.\n" % (carp_mesh))
        sys.exit(-1)

    # convert
    carp2Gmsh(sys.argv[1], sys.argv[2]) 

# end of main
