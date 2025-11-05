#!/usr/bin/env python

import sys, time, os
import numpy as np

def abaqusGetNodeset (inpMesh):  
    """
    Given a mesh in ABAQUS .inp file format, prints the nodeset section.
    ********************************** N O D E S E T S **************************
    *NSET, NSET=NS1
         101,     102,     103,     104,     105,     106,
         107,     108,     109,     110,     111,     112,
         113,     114,     115,     116,
    *NSET, NSET=NS2
         137,     138,     139,     140,     141,     142,
         143,     144,     145,     146,     147,     148,
         149,     150,     151,     152,
    **
    """
    
    inputFile = open(inpMesh)

    # find nodes section
    line = inputFile.readline()
    while(line.find("N O D E S E T S") == -1):
        line = inputFile.readline()

    line = inputFile.readline()
    num = 0
    nslst = []
    while(line.find("**") == -1):
        if (line.find("*NSET")!=-1):
            num = num + 1
            nslst.append([])
            line = inputFile.readline()
        ltemp = line.split(',')
        ltemp = map(int,ltemp[:-1])
        for ind in ltemp:
            nslst[num-1].append(ind)
        
        line = inputFile.readline()

    # print the nodeset's for my code
    for i,ns in enumerate(nslst):
        print("Nodeset %d, len %d" % (i,len(ns)))
        for node in ns:
            print('%d 0 0.0' % (node-1))
            print('%d 1 0.0' % (node-1))
            print('%d 2 0.0' % (node-1))

    # end

if __name__ == "__main__":

    if (len(sys.argv) < 2):
        print("\n Usage: abaqusGetNodeset <mesh.inp>\n")
        sys.exit(-1)

    # parse and check input
    inpfile = sys.argv[1]

    if (not os.path.isfile(inpfile)):
        print("\n Error: the input %s does not exist.\n" % (inpfile))
        sys.exit(-1)

    # convert
    abaqusGetNodeset(inpfile) 

# end of main
