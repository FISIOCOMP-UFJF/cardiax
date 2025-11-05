#!/usr/bin/env python

import os, sys, shutil
from sctools import carpElemLabelNumNode as getNumNode
from sctools import read_array_pts, read_array_elem, iseq, get_element_center

def fib2meshalyzer (filename):
    """
    Given a filename for a CARP mesh (.pts, .elem and .lon) this functions
    generates the .vec and .vpts files necessary to visualize vector data
    in meshalyzer.
    """
    ptsFile  = '%s.pts' % filename
    lonFile  = '%s.lon' % filename
    vecFile  = '%s.vec' % filename
    elemFile = '%s.elem' % filename
    vptsFile = '%s.vpts' % filename

    np_xyz   = read_array_pts (ptsFile)
    elemList = read_array_elem (elemFile)

    # copy filename.lon to filename.vec
    shutil.copy (lonFile, vecFile)   

    # create .vpts file
    f = open(vptsFile,'w')
    f.write('%d\n' % len(elemList))
    
    for e in elemList:
        numNode = getNumNode[e[0]]
        nodeList = [e[j+1] for j in xrange(numNode)]                        
        centroid = get_element_center(nodeList, np_xyz)
        f.write('%f %f %f\n' % (centroid[0], centroid[1], centroid[2]) )
    
    f.close()
    
    print ".vec and .vpts files were generated!"

# end of fib2meshalyzer

def main():
    
    args = sys.argv
    if len(args) < 2:
        print '\n Usage: python fib_meshalyzer <FILENAME>\n'
        exit(1)

    filename = args[1]    
    fib2meshalyzer(filename)
    
if __name__ == "__main__":
    
    #import profile
    #profile.run('main()')
    
    main()
