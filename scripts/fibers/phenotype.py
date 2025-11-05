#!/usr/bin/env python
# 
# This script is used to configure a CUBE mesh (.pts .elem) 
# for electromechanical simulations using fiber directions
# and also different cell phenotypes
#
# Bernardo M. Rocha
#
import sys, time, os
import numpy as np
from math import pi, sin, cos
from sctools import read_array_pts, read_array_elem

EPI, MCELL, ENDO = 0, 1, 2

def usage():
    print("""
    Usage: carp2gmsh <carp_mesh> <pEndo> <pM> <pEpi> <fibers>

    Options:
      carp_mesh: base name of the mesh
      pEndo: percentage of endocardial cells (default: 0.12)
      pM: percentage of M-cells              (default: 0.60)
      pEpi: percentage of epicardial cells   (default: 0.28)
      fibers: 1 to create transmural fiber directions (f0,s0,n0)
    """)

def createPhenotypes(xnodes, pendo=0.12, pmcell=0.6, pepi=0.28):
    """
    Given the points of a mesh, create 3 layers of the
    cardiac cell phenotypes: endo, M and epi and returns
    as an array.
    """
    n = len(xnodes)
    celltypes = np.zeros(n)

    xx = xnodes[:,0]
    x0 = xx.min()
    x1 = xx.max()
    xl = x1-x0

    for i in range(n):
        x = xnodes[i,0]
        xr = (x-x0)/xl
        if (xr < pendo):
            celltypes[i] = ENDO
        elif (xr >= pendo and xr <= pendo+pm):
            celltypes[i] = MCELL
        elif (xr >= pendo+pm and xr <= 1.0):
            celltypes[i] = EPI

    return celltypes        

def createFibers(xnodes, alfa=-60):
    """
    Create the fiber directions in the transmural setting
    with angle +alfa and -alfa at the endocardial and
    epicardial surfaces (x axis), respectively.
    """
    n = len(xnodes)   
    fvec, svec, nvec = np.zeros((n,3)), np.zeros((n,3)), np.zeros((n,3))

    # sheet directions
    for i in range(n):
        svec[i,0], svec[i,1], svec[i,2] = 1.0, 0.0, 0.0

    # find angle in rad
    alfa = (alfa*pi)/180.0
    alpi = 2*alfa 
    print("Endocardial %f" % alfa)
    print("Epicardial %f" % alpi)

    # fiber directions
    x0 = xnodes[:,0].min()
    x1 = xnodes[:,0].max()
    xl = x1-x0
    for i in range(n):
        x = xnodes[i,0]
        xr = (x-x0)/xl
        alfaNew = alfa - alpi*xr
	print("x=%f alfa=%f" % (x,(alfaNew*180)/pi)) 
        f0 = np.array([0.0, sin(alfaNew), cos(alfaNew)])
        f0 = f0/np.linalg.norm(f0)
        fvec[i,:] = f0[:]

    # normal directions
    for i in range(n):
        n0 = np.cross(fvec[i],svec[i])
        n0 = n0/np.linalg.norm(n0)
        nvec[i,:] = n0[:]

    return fvec, svec, nvec

# Main program
if __name__ == "__main__":

    # default values
    makefibers = False
    pendo, pm, pepi = 0.12, 0.60, 0.28
    
  # input
    if (len(sys.argv) < 2):
        usage()
        sys.exit(-1)

    # parse and check input
    carpMesh = sys.argv[1]
    
    if (len(sys.argv)==5):
        pendo = int(sys.argv[2])/100.
        pm = int(sys.argv[3])/100.
        pepi = int(sys.argv[4])/100.

    if (len(sys.argv)==6):
        if(int(sys.argv[5])==1):
            makefibers = True
    
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

    # Calculate phenotypes
    celltypes = createPhenotypes(pts, pendo, pm, pepi)

    # Write phenotypes to file
    typname = carpMesh + '.typ'
    f = open(typname,'w')
    f.write('%d\n' % numPts)
    for i in range(numPts):
        f.write('%d\n' % celltypes[i])
    f.close()

    # Are we going to make the fibers?
    if makefibers:
        # Compute fiber directions
        f0,s0,n0 = createFibers(pts)

        # Write fibers to files
        fibname = carpMesh + '.fib'
        print("Creating %s" % fibname)
        f = open(fibname,'w')
        f.write('0\n')
        for i in range(numPts):
            f.write(' %f %f %f' % (f0[i,0], f0[i,1], f0[i,2]))
            f.write(' %f %f %f' % (s0[i,0], s0[i,1], s0[i,2]))
            f.write(' %f %f %f' % (n0[i,0], n0[i,1], n0[i,2]))
            f.write('\n')
        f.close()

