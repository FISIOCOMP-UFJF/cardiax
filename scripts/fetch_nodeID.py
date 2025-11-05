#!/usr/bin/env python

import sys, pdb
import numpy as np
from sctools import read_array_pts

def get_node_ID (pts_file, x=0.0, y=0.0, z=0.0):
    """
    Given x,y and z coordinates return the ID of the node if it exists,
    if it does not exists returns -1.
    """
    xyz = read_array_pts(pts_file)

    for i in xrange(len(xyz[:,0])):
        if ((xyz[i][0] == x) and (xyz[i][1] == y) and (xyz[i][2] == z)):
            #print i, xyz[i]
            return i

    return -1

# end of get_node_ID
    
if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print "\n Usage: get_node_ID <pts_file> <x> <y> <z>\n"
        sys.exit(-1)

    pts_file = sys.argv[1]
    xcoord = float(sys.argv[2])
    ycoord = float(sys.argv[3])
    zcoord = float(sys.argv[4])

    index = get_node_ID (pts_file, xcoord, ycoord, zcoord)

    if index != -1:
        print ' Node ID: %d' % (index)
    else:
        print ' Node not found'

# end of file
    

