#!/usr/bin/env python

import os, sys, subprocess, pdb
import numpy as np
#from scipy.io import write_array
#from scipy.io import read_array
#from scipy.io.numpyio import fread
from scipy import array as sarray
from subprocess import Popen,PIPE

"""
A set of useful functions for reading text and binary files, executing command 
lines, generate integer and float sequences and so on.

Bernardo M. Rocha, 2008
"""

# interface to CARP elements
carp_map_label_node= {'Tt': 4, 'Hx': 8, 'Oc': 6, 'Py': 5,
                      'Pr': 6, 'Qd': 4, 'Tr': 3, 'Ln': 2}

def iseq (start=0, stop=None, inc=1):
    """
    Purpose: Generate integer from start to (and including) stop
    with increment of inc. Alternative to range/xrange.
    """
    if stop == None:
        stop = start; start = 0; inc = 1
    return xrange(start, stop+inc, inc)

def sequence (start, stop, inc):
    """
    Purpose: Generate float sequence from start to (and including) stop
    with increment of inc. Alternative to arange.
    """
    return np.arange(start, stop+inc, inc)

def read_array_pts (ptsFile):
    """
    Purpose: Function to read a .pts file from CARP simulator, where the first
    line contains the number of nodes and the remaining lines contain
    the array of nodes
    
    Example:
    from sctools import read_array_pts
    nn, narray = read_array_pts("t0010um_i.pts")
    print nn
    print narray[:,0]
    print narray[:,1]
    print narray[:,2]
    """
    try:
        f = open(ptsFile)
        lines = f.readlines()
        f.close()
    except IOError:
        print("IOError: File %s not found." % ptsFile)
        sys.exit(-1)
    
    # number of nodes
    numNodes = int(lines[0]); del lines[0]
    nodesTmp = "mytemp.txt"

    try:
        fnew = open(nodesTmp,"w")
        for i in lines:
            fnew.write(i)
        fnew.close()
    except IOError:
        print("IOError: File %s not found." % ptsFile)
        exit(-1)

    # numpy array of nodes
    nodes = np.loadtxt(nodesTmp)
    
    if(numNodes != len(nodes)):
        print(" error read_array_pts(): the size of the array doesn't match")
        exit(-2)
    
    os.remove(nodesTmp)
            
    return nodes.copy()

# end of read_array_pts

def read_array_elem (elemFile):
    """
    Purpose: Function to read a .elem file from CARP, where the first
    line contains the number of elements and the remaining lines contain
    the element conectivity
    
    Returns a list
    """
    f = open(elemFile)
    
    # extract number of elements
    header = f.readline()
    numElements = int(header.strip())
    
    # extract element list
    elemList = [(line.strip()).split(' ') for line in f.readlines()]
    
    f.close()
    
    return elemList

# end of read_array_elem

def read_binary_array (binfile, read_type):
    """
    Open a binary file read an array of read_type datatype
    and return the numpy an array
    PS: read_type is a character in 'cb1silfdFD' (PyArray types)
    """
    file = open(binfile,mode='rb')
    size = fread(file,1,'i')
    print(size)
    data = fread(file,size,read_type)
    file.close()
    return data

# end of read_binary_array

def read_map_array (mapfile):
    """
    """
    f = open(mapfile)
    size = int( f.readline().strip() )   
    maplist = [(line.strip()).split(' ') for line in f.readlines()]
    f.close()
    
    return np.array(maplist,dtype=np.int)

# end of read_map_array

def run_command_line (command, output=False):
    """
    Purpose: Run a command line and capture the stdout and stderr
    """   
    proc = subprocess.Popen(command,shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
    stdout_value, stderr_value = proc.communicate()	
    #if output:
	#   print('\tcombined output:', repr(stdout_value))

# end of run_command_line

def run_command_line_old (command, output=False):
    """
    Purpose: Run a command line and capture the stdout and stderr
    Deprecated: the use of popen2.popen4 will not be supported in Python 3
    """
    r,w = popen2.popen4(command)
    out = r.readlines()
    if output: print(out)
    r.close()
    w.close()

# end of run_command_line

def check_path (prog_name):
    """
    Purpose: find if the user has the prog_name in his $PATH
    """
    prog = False    
    pathList = os.environ.get('PATH').split(':')
    for path in pathList:
        if not os.access(path,os.X_OK):
            continue
        binaryList = os.listdir(path)
        binaryList.sort()
        if prog_name in binaryList:
            prog = True
            #if DEBUG: print "%s path=%s" % (prog_name,
            #os.path.join(path,carpBinary))
    if not prog:
        print("Error: %s was not found in your $PATH" % (prog_name))
        exit(-1)
        
# end of check_path

def get_element_center (nodeList, xyz):
    """
    Get the centroid of an element given the local node list and the array of
    global coordinates
    """
    x = np.zeros(len(nodeList))
    y = np.zeros(len(nodeList))
    z = np.zeros(len(nodeList))
    
    for i in xrange(len(nodeList)):
        index = nodeList[i]
        x[i] = xyz[index,0]
        y[i] = xyz[index,1]
        z[i] = xyz[index,2]
    
    return array([np.mean(x),np.mean(y),np.mean(z)])

# end of getElementCenter

def intersect_array1d (a,b,rows):
    """
    Finds set intersection of two vectors.    
    Returns: c the intersect vector and index vectors ia and ib such
    that c = a(ia) and c = b(ib).
    """    
    c  = np.intersect1d(a,b)
    ma = np.setmember1d(a,b)
    mb = np.setmember1d(b,a)
    ia = np.nonzero(ma)[0]
    ib = np.nonzero(mb)[0]
    return c, ia, ib

# end of intersect_array

def intersect_nodes (fpts1, fpts2, fout):
    """
    Given two .pts files, this function computes (maps) the common nodes
    to both files. The output is a file (fout) with the matching nodes indexes.
    """
    a = read_array_pts(fpts1)
    b = read_array_pts(fpts2)

    ia = np.logical_or.reduce(np.logical_and.reduce(a == b[:,None], axis=2))
    ib = np.logical_or.reduce(np.logical_and.reduce(b == a[:,None], axis=2))
    index_ia = np.array(np.nonzero(ia), dtype=np.int).squeeze()
    index_ib = np.array(np.nonzero(ib), dtype=np.int).squeeze()
    
    f = open(fout, 'w')
    f.write('%d\n' % (np.shape(index_ia)[0]))
    for i in xrange(len(index_ia)):
        #print index_ia[i], index_ib[i]
        f.write('%d %d\n' % (index_ia[i], index_ib[i]))
    f.close()
    
# end of intersect_nodes

if __name__ == "__main__" :
	a = '/data/sim/macrofem/papillary75um/hyb/papillary_75um.pts'
	b = '/data/sim/macrofem/papillary75um/tet/papillary_75um_tt.pts'
	o = '/home/rocha/mapa.txt'
	intersect_nodes(a,b,o)
