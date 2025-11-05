#!/usr/bin/env python
#
# Script to convert a mesh from Gmsh to CARP format
# Bernardo M. Rocha, 2014
#
import sys, time, os
import numpy as np

CHANGE_BBOX = False

def bbox(coords):
    xmin,xmax = coords[:,0].min() , coords[:,0].max()
    ymin,ymax = coords[:,1].min() , coords[:,1].max()
    zmin,zmax = coords[:,2].min() , coords[:,2].max()
    print('[%f , %f] x [%f , %f] x [%f, %f]' % (xmin,xmax,ymin,ymax,zmin,zmax))
    return (xmin,xmax,ymin,ymax,zmin,zmax)

def gmsh2carp (gmshMesh, outputMesh):
    
    ptsFile = outputMesh + '.pts'
    elemFile = outputMesh + '.elem'
    fibFile = outputMesh + '.fib'

    # read header from .msh file
    f = open(gmshMesh)
    line = f.readline() # $MeshFormat
    line = f.readline() # $2.1 0 8
    line = f.readline() # $EndMeshFormat
    
    # read nodes and write .pts file
    fpts = open(ptsFile,'w')
    line = f.readline() # $Node
    line = f.readline() # num_node
    num_nodes = int(line)
    fpts.write("%d\n" % (num_nodes))

    vpts = np.zeros((num_nodes,3))

    for i in range(num_nodes):
        line = f.readline().split(" ")
        node_id = int(line[0])
        x, y, z = float(line[1]), float(line[2]), float(line[3])
        vpts[i,0] = x
        vpts[i,1] = y
        vpts[i,2] = z

    if(CHANGE_BBOX):
        vb = bbox(vpts)    
        vpts[:,0] = vpts[:,0] + (-1) * vb[0]
        vpts[:,1] = vpts[:,1] + (-1) * vb[2]
        vpts[:,2] = vpts[:,2] + (-1) * vb[4]

    print('Bounding box')
    vb = bbox(vpts)

    for i in range(num_nodes):
        x,y,z = vpts[i,:]
        fpts.write("%f %f %f\n" % (x,y,z))

    line = f.readline() # $EndNodes
    fpts.close()
    print(" Reading nodes: Done.")
    
    # read elements and
    felem = open(elemFile,'w')
    line = f.readline() # $Elements
    line = f.readline() # num_elements

    num_elements = int(line)   
    elements = 0
    nodes = []

    for i in range(num_elements):
        line = f.readline().split(" ")
        elem_id = int(line[0])
        elem_type = int(line[1])

        if (elem_type == 1): # line segment
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[2+j]))            
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            nodes.append( (n1,n2) )
            #felem.write("Tr %d %d %d 1\n" % (n1,n2,n3))   
        elif (elem_type == 2): # triangle
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[3+j]))            
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            n3 = int(line[5+num_tags])
            nodes.append( (n1,n2,n3) )
            #felem.write("Tr %d %d %d 1\n" % (n1,n2,n3))   
        elif (elem_type == 5): # hexahedra
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[3+j]))
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            n3 = int(line[5+num_tags])
            n4 = int(line[6+num_tags])
            n5 = int(line[7+num_tags])
            n6 = int(line[8+num_tags])
            n7 = int(line[9+num_tags])
            n8 = int(line[10+num_tags])
            nodes.append( (n1,n2,n3,n4,n5,n6,n7,n8) )
        elif (elem_type == 4): # tetrahedra
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[3+j]))
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            n3 = int(line[5+num_tags])
            n4 = int(line[6+num_tags])
            t = (n1,n2,n3,n4)
            nodes.append( t )

    line = f.readline() # $EndElements
    f.close()
    print(" Reading elements: Done.")

    print(" Number of nodes: %d" % num_nodes)
    print(" Number of elements: %d" % num_elements)

    #elem_type = 4
    #nelem = 262145

    felem.write("%d\n" % (elements))    
    for i in range(elements):
        if(elem_type == 1):
            n1,n2 = nodes[i]
            felem.write("Ln %d %d 1\n" % (n1-1,n2-1))
        elif(elem_type == 2):
            n1,n2,n3 = nodes[i]
            felem.write("Tr %d %d %d 1\n" % (n1-1,n2-1,n3-1))
        elif(elem_type == 5):
            n1,n2,n3,n4,n5,n6,n7,n8 = np.array(nodes[i])-1
            felem.write("Hx %d %d %d %d %d %d %d %d 1\n" % (n1,n2,n3,n4,n5,n6,n7,n8))
        elif(elem_type == 4):
            n1,n2,n3,n4 = np.array(nodes[i])-1
            felem.write("Tt %d %d %d %d 1\n" % (n1,n2,n3,n4))

    felem.close()
    
    # create simple .fib file
    ffib = open(fibFile,'w')
    for i in range(elements):
        ffib.write("%f %f %f %f %f %f %f %f %f\n" % (1,0,0,0,1,0,0,0,1))
    ffib.close()    

# end of gmsh2carp  

if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: gmsh2carp <gmsh_mesh> <carp_output>\n"); sys.exit(-1)

    # parse and check input
    gmsh_mesh = sys.argv[1]

    if (not os.path.isfile(gmsh_mesh)):
       print("\n Error: the input gmsh %s does not exist.\n" % (gmsh_mesh))
       sys.exit(-1)

    # convert
    gmsh2carp(sys.argv[1], sys.argv[2]) 

# end of main
