#!/usr/bin/env python

import sys, time, os
import numpy as np
from sctools import read_array_pts, read_array_elem

def flagshyp2gmsh (flagMesh, outputMesh):

    # open output file
    print "Saving %s" % outputMesh
    o = open(outputMesh, 'w')
    
    # read header from .msh file
    f = open(flagMesh)
    tline = f.readline()      # title
    etype = f.readline()      # element type
    etype = etype.split()[0]
    nnode = int(f.readline()) # num nodes
    
    o.write("$MeshFormat\n")
    o.write("2.1 0 8\n")
    o.write("$EndMeshFormat\n")

    # parse and write nodes block
    o.write("$Nodes\n")
    o.write("%d\n" % nnode)
    for i in range(nnode):
        line = f.readline().split()
        if (etype == 'quad4'):
            nid, x,y = int(line[0]), float(line[2]), float(line[3])
            o.write("%d %f %f 0.0\n" % (nid,x,y))
    o.write("$EndNodes\n")

    line = f.readline()
    nelem = int(line)   

    # parse and write elements block
    o.write("$Elements\n")
    o.write("%d\n" % nelem)
    for i in range(nelem):
        line = f.readline().split()
        if (etype == 'quad4'):
            eid = int(line[0])
            a,b,c,d = int(line[2]),int(line[3]),int(line[4]),int(line[5])
            o.write("%d 3 3 9 8 0 %d %d %d %d\n" % (eid,a,b,c,d))
    o.write("$EndElements\n")            

    f.close()
    o.close()

if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print "\n Usage: flagshyp2gmsh <flagshyp_mesh> <gmsh>\n"; sys.exit(-1)

    # parse and check input
    flag_mesh = sys.argv[1]

    if (not os.path.isfile(flag_mesh)):
       print "\n Error: the input mesh %s does not exist.\n" % (flag_mesh)
       sys.exit(-1)

    # convert
    flagshyp2gmsh(sys.argv[1], sys.argv[2]) 

# end of main
