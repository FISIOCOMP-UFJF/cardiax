#!/usr/bin/env python

import sys
import numpy as np

def writeNodesFile(X,Y,Z,filename):

    f = open(filename,"w")

    sizex = np.size(X)
    sizey = np.size(Y)
    sizez = np.size(Z)

    if (sizex != sizey or sizex != sizez or sizey != sizez):
        print " writeNodesFile error: dimensions do not match.\n"
        exit(1)      

    size = sizex
    f.write("%d\n" % size)
    for i in xrange(size):
        f.write("%f %f %f \n" % (X[i],Y[i],Z[i]))

    f.close()
    
def writeElementsFile(filename):
    pass

def writeFibersFile(filename):
    pass

def writeVerticesFile(filename):
    pass


