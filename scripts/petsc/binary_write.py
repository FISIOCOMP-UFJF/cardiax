#!/usr/bin/env python

import sys, gzip, pdb, time
import numpy as np
import scipy.sparse as sparse
from scipy.io.numpyio import fwrite, fread
import matplotlib.pyplot as plt
from pyUtils.petsc_binary_read import petsc_binary_read

"""
Reads in PETSc binary file matrices or vectors
emits as scipy.sparse matrice in COO format or#
as a numpy vector

Description:
    1211216 to read a matrix
    1211214 to read a vector

Bernardo M. Rocha
"""

code_vector = 1211214
code_matrix = 1211216

def petsc_binary_write(A, filename):
    
    fd = open(filename,'wb')

    if np.rank(A) == 1:   # is vector
        m = np.shape(A)[0]
        n = 1
    elif np.rank(A) == 2: # is matrix
        m = np.shape(A)[0]
        n = np.shape(A)[1]

    if sparse.issparse(A):
        majic = 1.2345678910e-30
        diag = sparse.extract_diagonal(A)
        for i  in xrange(len(diag)):
            if diag[i] == 0:
                diag[i] = majic

        nz = sparse.spmatrix.getnnz(A)
        harr = np.array([1211216,m,n,nz])
        fwrite(fd, np.size(harr), harr, 'i')

        #nz_per_row = 

        #write nz_per_row

#    fwrite(fd,n_nz,'int32');  %nonzeros per row
#    [i,j,s] = find(A');
#    fwrite(fd,i-1,'int32');
#    for i=1:nz
#      if s(i) == majic
#        s(i) = 0;
#      end
#    end
#    fwrite(fd,s,'double');


    else:
        size = m * n
        harr = np.array([code_vector, size]) # header
        fwrite(fd, 2, harr, 'i', 1)
        fwrite(fd, size, A, 'd', 1)

    pdb.set_trace()

    fd.close()
       
# end of petsc_binary_write

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "\n Usage: petsc_binary_write [filename] \n"
        exit(-1)
    
    inputfile = str(sys.argv[1])

    out = petsc_binary_read(inputfile, 0)

    petsc_binary_write(out,'saida.bin')

    
