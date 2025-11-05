#!/usr/bin/env python
#
# Create and configure the .fib for the simple cube mesh
# Bernardo M. Rocha
#
import sys, time, os
import numpy as np
from sctools import read_array_pts, read_array_elem
from fiberUtil import rotateAroundAxis
from math import pi, ceil, cos, sin, tan

if __name__ == "__main__":

    if (len(sys.argv) < 2):
        print("\n Usage: configFiberCube <carp_mesh>\n")
        sys.exit(-1)

    # parse and check input
    carpMesh = sys.argv[1]

    if (not os.path.isfile(carpMesh+'.pts')) or (not os.path.isfile(carpMesh+'.elem')):
        print("\n Error: the input carpfile %s does not exist.\n" % (carpMesh))
        sys.exit(-1)

    ptsFile  = carpMesh + '.pts'
    elemFile = carpMesh + '.elem'

    pts  = read_array_pts (ptsFile)
    elem = read_array_elem (elemFile)

    numPts  = np.shape(pts)[0]
    numElem = np.shape(elem)[0]

    print('Number of points: %d' % numPts)
    print('Number of elements: %d' % numElem)

    # start
    x0 = pts[:,0].min()
    x1 = pts[:,0].max()
    xl = x1-x0

    nptx = ceil(pow(numPts,1./3.))
    npty, nptz = nptx, nptx
    sx = np.linspace(x0,x1,nptx)

    # angles for endo and epi
    aendo, aepi = 70.0, -70.0

    # define some vectors
    ax = np.array([1.,0.,0.]) 
    f0 = np.array([0.,0.,1.]) 

    # create fibers 
    fibfile = open('fibers.fib','w')
    fibfile.write('0\n')

    for i in range(numPts):
        x = pts[i,0]
        r = (x-x0)/xl        
        alfa = (2.0*pi*(1.0-r))/3.0 - (pi/4.0)
        f = rotateAroundAxis(f0, ax, alfa)
        f = f/np.linalg.norm(f)
        s = np.array([1.0, 0.0, 0.0])
        n = np.cross(f,s)
        n = n/np.linalg.norm(n)
        
        # transversely isotropic model
        #fibfile.write('%f %f %f\n' % (f[0],f[1],f[2]))
        
        # orthotropic model
        fibfile.write('%f %f %f '  % (f[0],f[1],f[2]))
        fibfile.write('%f %f %f '  % (s[0],s[1],s[2]))
        fibfile.write('%f %f %f\n' % (n[0],n[1],n[2]))
        
        ix = np.where(sx==x)
        
        

# end of main
