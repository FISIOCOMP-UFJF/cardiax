#!/usr/bin/env python
#
# Script to convert a mesh from Gmsh to Cardiax XML format
# Gmsh element types
# 1: 2-node line.
# 2: 3-node triangle.
# 3: 4-node quadrangle.
# 4: 4-node tetrahedron.
# 5: 8-node hexahedron.
# 6: 6-node prism.
# 7: 5-node pyramid.
#
# Bernardo M. Rocha, 2015
#
import sys, time, os
import numpy as np

gmsh_etype_dim = {'1': 1, '2': 2, '3': 2, '4': 3, '5': 3}

def gmsh2xml (gmshMesh, outputMesh, parFile=None):

    xmlfilename = outputMesh

    # read header from .msh file
    f = open(gmshMesh)
    line = f.readline() # $MeshFormat
    line = f.readline() # $2.1 0 8
    line = f.readline() # $EndMeshFormat

    # read nodes and write .pts file
    line = f.readline() # $Node
    line = f.readline() # num_node
    num_nodes = int(line)
    pts = np.zeros((num_nodes,3))
    for i in range(num_nodes):
        line = f.readline().split(" ")
        node_id = int(line[0])
        x, y, z = float(line[1]), float(line[2]), float(line[3])
        pts[i,0] = x
        pts[i,1] = y
        pts[i,2] = z

    line = f.readline() # $EndNodes
    print("Reading nodes...done")

    # read elements and
    line = f.readline() # $Elements
    line = f.readline() # num_elements

    num_elements = int(line)
    elements = 0
    elems = []
    elem_types = []
    elem_marker = []

    for i in range(num_elements):
        line = f.readline().split(" ")
        elem_id = int(line[0])
        elem_type = int(line[1])
        elem_types.append(elem_type)

        if (elem_type == 1): #line
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[3+j]))
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            elems.append( (n1,n2) )
            #print("marker %d" % tags[0] )
            #print(tags)
            elem_marker.append( tags[0] )
        elif (elem_type == 2): # triangle
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[3+j]))
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            n3 = int(line[5+num_tags])
            elems.append( (n1,n2,n3) )
            #print("marker %d" % tags[0])
            elem_marker.append( tags[0] )
            #felem.write("Tr %d %d %d 1\n" % (n1,n2,n3))
        elif (elem_type == 3): # quad
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[3+j]))
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            n3 = int(line[5+num_tags])
            n4 = int(line[6+num_tags])
            elems.append( (n1,n2,n3,n4) )
            #print("marker %d" % tags[0] )
            elem_marker.append( tags[0] )
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
            elems.append( (n1,n2,n3,n4,n5,n6,n7,n8) )
            elem_marker.append( tags[-2] )
        elif (elem_type == 4): # tetrahedra
            elements = elements + 1
            num_tags = int(line[2])
            tags = []
            for j in range(num_tags):
                tags.append(int(line[3+j]))
            #print(tags[-2])
            n1 = int(line[3+num_tags])
            n2 = int(line[4+num_tags])
            n3 = int(line[5+num_tags])
            n4 = int(line[6+num_tags])
            t = (n1,n2,n4,n3)
            elems.append( t )
            elem_marker.append( tags[0] )

    line = f.readline() # $EndElements
    f.close()

    # make some considerations
    # gmsh element type 1-line, 2-tri, 3-quad, 4-tet, 5-hex, ....
    elem_type = np.array(elem_types).max()
    if (elem_type <= 1 or elem_type > 5):
        print("Error: element type is not supported")
        sys.exit(1)

    # find out dimension
    num_dim = 1
    if (elem_type == 2 or elem_type == 3):   num_dim = 2
    elif (elem_type == 4 or elem_type == 5): num_dim = 3

    # read boundary elements section
    selem = []
    velem = []
    selem_markers = []

    for i in xrange(num_elements):
        etype = int(elem_types[i])
        eskey = str(etype)
        eldim = gmsh_etype_dim[eskey]
        marker = elem_marker[i]

        # then it is a boundary element
        if(eldim < num_dim):
            selem.append( np.array(elems[i][:],dtype=int) )
            selem_markers.append( elem_marker[i] )
        else:
            velem.append( np.array(elems[i][:],dtype=int) )

    num_elements = len(velem)
    num_boundary_elem = len(selem)

    print("Reading elements...done.")

    elem_xml = "line"
    if(elem_type == 2):
        elem_xml = "triangle"
    elif(elem_type == 3):
        elem_xml = "quadrilateral"
    elif(elem_type == 4):
        elem_xml = "tetrahedron"
    elif(elem_type == 5):
        elem_xml = "hexahedron"

    # print info
    print("Number of dimensions: %d" % num_dim)
    print("Number of nodes: %d" % num_nodes)
    print("Number of elements: %d" % num_elements)
    print("Number of boundary elements: %d" % num_boundary_elem)
    print("Element type: %s" % elem_xml)

    # reading parameters
    prescDispl = []
    problemtyp = None
    if parFile is not None:
        with open(parFile, "r") as fpar:
            searchlines = fpar.readlines()

        for i, line in enumerate(searchlines):

            if ("prescribed_displacement" in line and line[0:4] != "end_"):
                numPresc = int(searchlines[i+1])
                j = i+2
                for l in searchlines[j:]:
                    if ("end_prescribed_displacement" in l):
                        break
                    lm = l.split()
                    nodeId  = int(lm[0])
                    nodeDir = int(lm[1])
                    nodeVal = float(lm[2])
                    prescDispl.append( (nodeId, nodeDir, nodeVal) )

                if("problem_type" in line):
                    l = line.split()
                    problemtyp = l[1]

        fpar.close()
        print("Reading parameters file %s...done" % parFile)


    base_indices = []
    for i in range(len(selem)):
        if selem_markers[i] == 10:
            base_indices.append(selem[i][:] - 1)

    epi_indices = []
    for i in range(len(selem)):
        if selem_markers[i] == 40:
            epi_indices.append(selem[i][:] - 1)

    prescDispl = np.unique(base_indices)
    prescDisplEPI = np.unique(epi_indices)

    #
    # start writing the output file
    #

    # basic header
    outputFile = file(outputMesh, 'w')

#    if(len(prescDispl) > 0):
#        if(problemtyp is not None):
#            outputFile.write('<elasticity type="%s">\n' % problemtyp)
#        else:
#            outputFile.write('<elasticity>\n')
#        outputFile.write('  <prescribed_displacement>\n')
#        for i in range(len(prescDispl)):
#            node,d,v = prescDispl[i]
#            outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,d,v))
#        outputFile.write('  </prescribed_displacement>\n')
#        outputFile.write('</elasticity>\n')

    outputFile.write('<elasticity>\n')

    outputFile.write('  <prescribed_displacement>\n')
    for i in range(len(prescDispl)):
        #node,d,v = prescDispl[i]
        #outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,d,v))
           
        node = prescDispl[i]
        #outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,0,0))
        #outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,1,0))
        outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,2,0))

    for i in range(len(prescDisplEPI)):
        node = prescDisplEPI[i]
        outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,0,0))
        outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,1,0))
        outputFile.write('    <node id="%d" direction="%d" value="%f" />\n' % (node,2,0))


    outputFile.write('  </prescribed_displacement>\n')
    outputFile.write('</elasticity>\n')


    # write mesh
    outputFile.write('<mesh celltype="%s" dim="%d">\n' % (elem_xml, num_dim))

    # nodes section
    outputFile.write('  <nodes size="%d">\n' % (num_nodes))
    for i in xrange(num_nodes):
        outputFile.write('    <node id="%d" ' % (i) )
        outputFile.write('x="%f" y="%f" z="%f" />\n' % (pts[i,0],pts[i,1],pts[i,2]))
    outputFile.write('  </nodes>\n')

    # elements section
    outputFile.write('  <elements size="%d">\n' % (num_elements))
    for i in xrange(len(velem)):
        #etype = int(elem_types[i])
        #eskey = str(etype)
        #eldim = gmsh_etype_dim[eskey]
        # then it is a boundary element
        #if(eldim < num_dim):
        #    boundary.append( np.array(elems[i][:],dtype=int)-1 )
        outputFile.write('    <element id="%d" ' % (i))
        write_element(outputFile, velem[i][:] )
    outputFile.write('  </elements>\n')

    # boundary section
    # first find type of boundary element
    # then write to file
    if(elem_type == 3): belem_xml = "line"
    elif(elem_type == 4): belem_xml = "triangle"
    elif(elem_type == 5): belem_xml = "quadrilateral"

    outputFile.write('  <boundary celltype="%s" dim="%d">\n' % (belem_xml,num_dim-1))
    for i in range(len(selem)):
        outputFile.write('    <element id="%d" marker="%d" ' % (i,selem_markers[i]))
        write_element(outputFile, selem[i][:] )
    outputFile.write('  </boundary>\n')

    outputFile.write('</mesh>\n')
    outputFile.close()
    print("Done")

# ------------------------------------------------------------------------------

def write_element(out, conec):
    num_nodes = len(conec)
    for i in range(num_nodes):
        out.write('v%d="%d" ' % (i,conec[i]-1))
    out.write('/>\n')

# ------------------------------------------------------------------------------

if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: gmsh2xml <gmsh_mesh> <output_xml> [parameters.par]\n")
        sys.exit(-1)

    # parse and check input
    gmsh_mesh = sys.argv[1]

    if (not os.path.isfile(gmsh_mesh)):
       print("\n Error: the input gmsh %s does not exist.\n" % (gmsh_mesh))
       sys.exit(-1)

    par = None
    if(len(sys.argv) == 4):
        par = sys.argv[3]

    # convert
    gmsh2xml(sys.argv[1], sys.argv[2], par)

# end of main
