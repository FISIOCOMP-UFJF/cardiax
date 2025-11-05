#!/usr/bin/env python
#
# Script to convert Abaqus mesh to Gmsh (.inp ---> .msh)
# Bernardo M. Rocha, 2014
#

import sys, time, os
import numpy as np

def nodesSection(gmshfile, vpts):
    numPts = np.shape(vpts)[0]
    gmshfile.write('$Nodes\n')
    gmshfile.write('%d\n' % numPts)
    for i in range(numPts):
        if( np.shape(vpts)[1] == 2):
            gmshfile.write('%d %f %f 0.0\n' % (i+1,vpts[i,0],vpts[i,1]))
        elif ( np.shape(vpts)[1] == 3):
            gmshfile.write('%d %f %f %f\n' % (i+1,vpts[i,0],vpts[i,1],vpts[i,2]))
    gmshfile.write('$EndNodes\n')

def unique(seq):
    seen = set()
    for item in seq:
        if item not in seen:
            seen.add( item )
            yield item

def getElementNumNodes(elem, eltype):
    if  (elem=="STRI3"):
        return 3
    elif(elem=="C3D4" or elem=="S4R"):
        return 4
    elif(elem=="C3D8R"):
        return 8
    #else:
    #    return int(eltype[k][3])

def abaqus2Gmsh (inpMesh, outputMesh):
    """
    Given a mesh in ABAQUS .inp file format, convert it to Gmsh's format.
    It creates only nodes and elements sections.
    """
    ptslist, elelist =  list(), list()
    etagdict  = {} # holds { (elementId, surfId), ... }
    etagcount = 0
    eltype = []
    numElem = []

    print("Parsing input file")

    with open(inpMesh, "r") as f:
        searchlines = f.readlines()

    for i, line in enumerate(searchlines):

        if ("*NODE" in line):
            for l in searchlines[i+1:]:
                if ("**" in l): break
                lm = map(float,l.split(',')[1:])
                ptslist.append(lm)

        if ("*ELEMENT" in line):
            tmpelemlist = []
            elinfo = line.split(",")
            print("Found element type %s" % elinfo[1].strip()[5:])
            eltype.append( elinfo[1].strip()[5:] )
            for l in searchlines[i+1:]:
                if ("**" in l): break
                if ("*ELEMENT" in l): break
                enodes = map(int,l.split(',')[1:])
                elelist.append( enodes ) 
                tmpelemlist.append( enodes )
            numElem.append( len(tmpelemlist) )

        if ("ELSET, ELSET=" in line):
            print("Converting element set")
            etagcount = etagcount + 1
            elsetinfo = line.split("=")[1]
            #print(searchlines[i])
            #print(searchlines[i+1])
            for l in searchlines[i+1:]:
                if ("**" in l): break
                if ("*SURFACE" in l or "SS" in l):
                    pass
                else:
                    l = l.strip()
                    l = l.rstrip(',')
                    etag = filter(int,l.split(','))
                    etag = map(int,l.split(','))
                    for i in range(len(etag)):
                        #etaglist[etagcount-1].append(etag[i])
                        elemId = etag[i]
                        etagdict[elemId] = etagcount

    f.close()
    print('Parsing Done\n')

    # split element list in case we have surface elements
    nvElem = 0
    nsElem = 0

    nl0 = len(elelist[0])
    idxsplit = -1
    for i in range(len(elelist)):
        nl = len(elelist[i])
        if(nl != nl0):
            idxsplit = i
            break

    element_list = elelist[idxsplit]

    print("Index to split element list: %d" % (idxsplit))

    if(idxsplit != -1):
        slst = elelist[:idxsplit]
        vlst = elelist[idxsplit:]
    else:
        print("No surface elements found.")
        slst = []
        vlst = elelist

    # create numpy arrays
    print('Creating arrays')
    vpts  = np.array(ptslist)
    velem = np.array(vlst, dtype=np.int32)
    selem = np.array(slst, dtype=np.int32)

    # write output file
    outputFile = open(outputMesh, 'w')
    outputFile.write('$MeshFormat\n2.0 0 8\n$EndMeshFormat\n')
    nodesSection(outputFile, vpts)

    # elements section
    # first volume elements, then surface elements
    totalElem = sum(numElem)
    #totalElem = numElem[0]
    #if(len(numElem)>1): 
    #    totalElem += numElem[1]

    print("Number of nodes: %d" % len(vpts))
    print("Number of volume elements: %d" % len(vlst))
    print("Number of surface elements: %d" % len(slst))
    print("Total number of elements: %d" % totalElem)

    # pre-process eltype list to remove duplicates
    #eltype = list( unique(eltype) )

    outputFile.write('$Elements\n')
    outputFile.write('%d\n' % totalElem)

    nsurfelem = len(slst)
    hasSurfel = False
    if nsurfelem > 0: hasSurfel = True

    # save VOLUME first then SURFACE elements
    # invert lists if necessary
    if(hasSurfel):
        eltype.reverse()
        numElem.reverse()
    print("Element types: %s" % eltype)
    print("Number of elements of each type: %s" % numElem)

    # volume elements AND surface elements
    if (hasSurfel):
        if(len(eltype) > 0):
            for k in range(len(eltype)):
                elem = eltype[k]

                if  (elem=="STRI3"):
                     numNodesElem=3
                elif(elem=="C3D4" or elem=="S4R"):
                    numNodesElem=4
                elif(elem=="C3D8R"):
                    numNodesElem=8
                else:
                    numNodesElem=int(eltype[k][3])

                print(" Element type = %s" % elem)
                print("  Number of elements = %d" % numElem[k])
                print("  Number of nodes per element = %d" % numNodesElem)

                for i in range(numElem[k]):
                    surfindex = 10 + k
                    if (i+1 in etagdict.keys()):
                        surfindex = etagdict[i+1]*10

                    e = i
                    if(k==1):   e = i + numElem[0]
                    elif(k==2): e = i + numElem[0] + numElem[1]

                    if (numNodesElem == 8):
                        outputFile.write('%d 5 2 %d 2 ' % (e+1,surfindex))
                    elif (numNodesElem == 4 and elem=='S4R'):
                        outputFile.write('%d 3 2 %d 2 ' % (e+1,surfindex))
                    elif (numNodesElem == 4 and elem=='C3D4'):
                        outputFile.write('%d 4 2 %d 2 ' % (e+1,surfindex))
                    elif (numNodesElem == 3):
                        outputFile.write('%d 2 2 %d 2 ' % (e+1,surfindex))

                    if (k == 0):   # volume elements first
                        outputFile.write('%s\n' % ' '.join(map(str,vlst[i])))
                    elif (k == 1): # then surface elements
                        outputFile.write('%s\n' % ' '.join(map(str,slst[i])))

        # end elemtype loop

        #outputFile.write('$EndElements\n')
        #outputFile.close()
        #print("Done")

    # volume elements only - NO SURFACE elements
    else:
        print("No surface elements")
        if(len(eltype) > 0):
            for k in range(len(eltype)):
                elem = eltype[k]

                numNodesElem = getElementNumNodes(elem, eltype)

                print(" Element type = %s" % elem)
                print("  Number of elements = %d" % numElem[k])
                print("  Number of nodes per element = %d" % numNodesElem)

                for i in range(numElem[k]):
                    volindex = 10
                    e = i

                    if (numNodesElem == 8):
                        outputFile.write('%d 5 2 %d 2 ' % (e+1,volindex))
                    elif (numNodesElem == 4 and elem=='S4R'):
                        outputFile.write('%d 3 2 %d 2 ' % (e+1,volindex))
                    elif (numNodesElem == 4 and elem=='C3D4'):
                        outputFile.write('%d 4 2 %d 2 ' % (e+1,volindex))
                    elif (numNodesElem == 3):
                        outputFile.write('%d 2 2 %d 2 ' % (e+1,volindex))

                    outputFile.write('%s\n' % ' '.join(map(str,vlst[i])))
        # end elemtype loop
    # end of

    outputFile.write('$EndElements\n')
    outputFile.close()
    print("Done")

# end of abaqus2Gmsh

if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: abaqus2gmsh <mesh.inp> <mesh.msh>\n")
        sys.exit(-1)

    # parse and check input
    inpfile = sys.argv[1]

    if (not os.path.isfile(inpfile)):
        print("\n Error: the input %s does not exist.\n" % (inpfile))
        sys.exit(-1)

    # convert
    abaqus2Gmsh(inpfile, sys.argv[2]) 

# end of main
