#!/usr/bin/env python

import sys, gzip, pdb, time
import numpy as np
from scipy import sparse
from scipy.io.numpyio import fwrite, fread
import matplotlib.pyplot as plt

"""
Reads in PETSc binary file matrices or vectors
emits as scipy.sparse matrice in COO format or#
as a numpy vector

Description:
    1211216 to read a matrix
    1211214 to read a vector

Bernardo M. Rocha
"""

def petsc_binary_read(filename, spymatrix):
    
    fd = open(filename,'rb');
    datatype = 'i'
    header = fread(fd, 2, datatype, datatype, 1)

    if np.size(header) <= 0:
        print 'File does not have that many items'
        sys.exit(1)

    if header[0] == 1211216:  # read matrix
        m      = header[1]
        header = fread(fd, 2, datatype, datatype, 1)
        n      = header[0]
        nz     = header[1]
        nnz    = fread(fd, m, datatype, datatype, 1)

        j   = fread(fd, nz, 'i', 'i', 1) + 1 # columns
        s   = fread(fd, nz, 'd', 'd', 1)     # data
        i   = np.ones((nz))                  # pointers

        cnt = 0
        for k in xrange(m):
            next          = cnt + nnz[k] - 1
            i[cnt:next+1] = k * np.ones((nnz[k]))
            cnt           = next + 1

        # adjust to python-zero based indexing
        j = j - 1

        A = sparse.coo_matrix( (s,(i,j)), (m,n) )

        if spymatrix == 1:
            plt.spy(A,marker='.')
            plt.show()

        fd.close()
        return A

    elif header[0] == 1211214: # read vector
        m = header[1]
        datatype = 'd'
        v = fread(fd, m, datatype, datatype, 1)
        
        fd.close()
        return v
       
# end of petsc_binary_read

if __name__ == "__main__":

    spymatrix = 0

    if len(sys.argv) < 2:
        print "\n Usage: petsc_binary_read [filename] [plot] \n"
        exit(-1)
    
    inputfile = str(sys.argv[1])
    if len(sys.argv) == 3:
        spymatrix = int(sys.argv[2])
    
    out = petsc_binary_read(inputfile, spymatrix)

    print out
    
