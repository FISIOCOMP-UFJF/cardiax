#!/usr/bin/env python

import sys, time, os, os.path
import numpy as np
from sctools import read_array_pts, read_array_elem, read_elem_type

def carp2xml (carpMesh, outputMesh):

    ptsFile  = carpMesh + '.pts'
    elemFile = carpMesh + '.elem'
    fibFile  = carpMesh + '.fib'
    stimFile = carpMesh + '.stim'
    
    # find element type
    elem_name = read_elem_type(elemFile)
    pts  = read_array_pts (ptsFile)
    elem = read_array_elem (elemFile)

    elem_xml = None
    num_dim  = None
    num_pts  = np.shape(pts)[0]
    num_elem = np.shape(elem)[0]

    # find out the corresponding XML element type
    if(elem_name == "Ln"): 
        num_dim, elem_xml = 1, "line"
    elif(elem_name == "Tr"):
        num_dim, elem_xml = 2, "triangle"
    elif(elem_name == "Qd"):
        num_dim, elem_xml = 2, "quadrilateral"
    elif(elem_name == "Tt"):
        num_dim, elem_xml = 3, "tetrahedron"
    elif(elem_name == "Hx"):
        num_dim, elem_xml = 3, "hexahedron"

    # print info
    print("Number of nodes: %d" % num_pts)
    print("Number of elements: %d" % num_elem)
    print("Element type: %s" % elem_xml)

    # basic header
    outputFile = file(outputMesh, 'w')
    outputFile.write('<?xml version=\"1.0\"?>\n')    
    outputFile.write('<mesh celltype="%s" dim="%d">\n' % (elem_xml, num_dim))

    # nodes section
    outputFile.write('  <nodes size="%d">\n' % (num_pts))
    for i in xrange(num_pts):
        outputFile.write('    <node id="%d" ' % (i) ) 
        outputFile.write('x="%f" y="%f" z="%f" />\n' % (pts[i,0],pts[i,1],pts[i,2]))
    outputFile.write('  </nodes>\n')

    # elements section
    outputFile.write('  <elements size="%d">\n' % (num_elem))
    for i in xrange(num_elem):
        etype = str(elem[i][0])
        outputFile.write('    <element id="%d" ' % (i))
        write_element(outputFile, np.array(elem[i][1:-1],dtype=int) )
    outputFile.write('  </elements>\n')

    # reading fibers
    if( os.path.isfile(fibFile) ):
        print("Reading fibers from %s" % fibFile)
        fvec = open(fibFile,"r")
        fibtype = int(fvec.readline())
        fibname = None
        fvec.close()    
        vecs = np.loadtxt(fibFile, skiprows=1)
        fibsize = np.shape(vecs)[0]
        if(fibsize != num_elem):
            print("Error: number of fiber vectors must match the number of elements")
            sys.exit(1)
        if(fibtype==1):
            fibname = "fiber_transversely_isotropic"
        else:
            fibname = "fiber_orthotropic"
        print("Fiber model: %s" % fibname)
        print("Fiber vectors: %d" % fibsize)

        # element_data section
        outputFile.write('  <element_data type="%s">\n' % (fibname))
        for i in xrange(fibsize):
            outputFile.write('    <element id="%d">\n' % (i))
            write_vecs(outputFile, vecs[i,:])
            outputFile.write('    </element>\n')
        outputFile.write('  </element_data>\n')        
    outputFile.write('</mesh>\n')


    # electrophysiology parameters
    outputFile.write('<electrophysiology>\n')
    if( os.path.isfile(stimFile) ):
        sfile = open(stimFile, "r")
        nstim = int(sfile.readline())
        vstim = []
        for i in range(nstim):
            l = sfile.readline()
            vstim.append(l)
        sfile.close()
        outputFile.write('  <stimuli number=\"%d\">\n' % nstim)
        for i in range(nstim):
            write_stimuli(outputFile, vstim[i])
        outputFile.write('  </stimuli>\n')
    outputFile.write('</electrophysiology>\n')
    
    outputFile.close()
    print("Done")

# ------------------------------------------------------------------------------    

def write_element(out, conec):
    num_nodes = len(conec)
    for i in range(num_nodes):
        out.write('v%d="%d" ' % (i,conec[i]))
    out.write(' />\n')

# ------------------------------------------------------------------------------

def write_vecs(out, vec):
    size = np.shape(vec)[0]
    # transversely isotropic - fiber only
    if(size == 3):
        out.write('        <fiber>%f,%f,%f</fiber>\n' % (vec[0],vec[1],vec[2]))
    elif(size == 9):
        out.write('        <fiber>%f,%f,%f</fiber>\n' % (vec[0],vec[1],vec[2]))
        out.write('        <sheet>%f,%f,%f</sheet>\n' % (vec[3],vec[4],vec[5]))
        out.write('        <normal>%f,%f,%f</normal>\n' % (vec[6],vec[7],vec[8]))
        # TERMINAR DE IMPLEMENTAR
        pass

# ------------------------------------------------------------------------------

def write_stimuli(out, s):
    stim = s.split()
    s  = float(stim[0])
    d  = float(stim[1])
    v  = float(stim[2])
    x0 = float(stim[3])
    x1 = float(stim[4])
    y0 = float(stim[5])
    y1 = float(stim[6])
    z0 = float(stim[7])
    z1 = float(stim[8])
    out.write('      ')
    out.write('<stim')
    out.write(' start=\"%.2f\" duration=\"%.2f\" value=\"%.4f\" ' % (s,d,v))
    out.write('x0=\"%.2f\" x1=\"%.2f\" ' % (x0,x1))
    out.write('y0=\"%.2f\" y1=\"%.2f\" ' % (y0,y1))
    out.write('z0=\"%.2f\" z1=\"%.2f\" ' % (z0,z1))
    out.write('/>\n')
    

# ------------------------------------------------------------------------------
    
if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: carp2xml <carp_mesh> <xml_output>\n")
        sys.exit(-1)

    # parse and check input
    carp_mesh = sys.argv[1]

    if (not os.path.isfile(carp_mesh+'.pts')) or (not os.path.isfile(carp_mesh+'.elem')):
        print("\n Error: the input carpfile %s does not exist.\n" % (carp_mesh))
        sys.exit(-1)

    # convert
    carp2xml(sys.argv[1], sys.argv[2]) 

# end of main
# ------------------------------------------------------------------------------
