#!/usr/bin/env python

import sys, time, os
import numpy as np

def febio2Gmsh (febMesh, outputMesh):

    f = open(febMesh)
    line = f.readline()
    while(line.find("<Geometry>") == -1):
        line = f.readline()
       
    # read points
    lpts = []
    line = f.readline() # <Nodes>
    line = f.readline()
    while(line.find("</Nodes>") == -1):
        temp = line.split(',')       
        xt = temp[0].split()

        if(len(xt) == 3):
            x = float(xt[2])
        else:
            xt = xt[1].split('>')
            x = float(xt[1])

        y = float(temp[1])
        zi = temp[2].index("</node>")
        zt = temp[2][:zi]
        z = float(zt)       
        lpts.append( [x,y,z] )
        line = f.readline()

    num_pts = len(lpts)
    vpts = np.array(lpts)

    # read elements
    lelem = []
    line = f.readline() # </Nodes>
    line = f.readline() # <Elements>
    while(line.find("</Elements>") == -1):
        temp = line.split(',')
        eltag  = temp[0].split()
        eltype = eltag[0][1:]
        
        if (eltype == 'hex8'):            
            node1 = int(eltag[-1])
            node2 = int(temp[1])
            node3 = int(temp[2])
            node4 = int(temp[3])
            node5 = int(temp[4])
            node6 = int(temp[5])
            node7 = int(temp[6])
            ti = temp[7].find("</hex8>")
            node8 = int(temp[7][:ti])
            lelem.append([node1,node2,node3,node4,node5,node6,node7,node8])

        elif (eltype == 'tet4'):
            node1 = int(eltag[-1])
            node2 = int(temp[1])
            node3 = int(temp[2])
            ti = temp[3].find("</tet4>")
            node4 = int(temp[3][:ti])
            lelem.append([node1,node2,node3,node4])

        
        line = f.readline()

    num_elem = len(lelem)
    velem = np.array(lelem)
    
    # basic header
    outputFile = open(outputMesh, 'w')
    outputFile.write('$MeshFormat\n2.0 0 8\n$EndMeshFormat\n')

    # nodes section
    outputFile.write('$Nodes\n')
    outputFile.write('%d\n' % num_pts)
    for i in xrange(num_pts):
        outputFile.write('%d %f %f %f\n' % (i+1,vpts[i,0],vpts[i,1],vpts[i,2]))
    outputFile.write('$EndNodes\n')

    # elements section
    outputFile.write('$Elements\n')
    outputFile.write('%d\n' % num_elem)
    
    for i in xrange(num_elem):
        if (eltype == 'hex8'):
            outputFile.write('%d 5 2 99 2 ' % (i+1))
        elif (eltype == 'tet4'):
            outputFile.write('%d 4 2 99 2 ' % (i+1))
        outputFile.write('%s\n' % ' '.join(map(str,velem[i])))
        
    outputFile.write('$EndElements\n')
    outputFile.close()


if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: febio2gmsh <file.feb> <file.msh>\n")
        sys.exit(-1)

    # parse and check input
    feb_mesh = sys.argv[1]
    out_mesh = sys.argv[2]

    if (not os.path.isfile(feb_mesh)):
        print("\n Error: the input file %s does not exist.\n" % (feb_mesh))
        sys.exit(-1)

    # convert
    febio2Gmsh(feb_mesh, out_mesh)

# end of main
