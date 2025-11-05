#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
from scipy import sparse
from pyUtils.sctools import read_binary_array

#file_Ax = '/home/rocha/Desktop/numerik_und_simulation/FEM2D_60/sparse_vals.bin'
#file_Aj = '/home/rocha/Desktop/numerik_und_simulation/FEM2D_60/sparse_cols.bin'
#file_Ap = '/home/rocha/Desktop/numerik_und_simulation/FEM2D_60/sparse_ptrs.bin'

def usage():
    print "\n Usage: sparseCSRspy <vals_array> <cols_array> <ptrs_array>\n"

# end of usage

def plot_csr_matrix(file_Ax, file_Aj, file_Ap):

    Ax = read_binary_array (file_Ax,'f')
    Aj = read_binary_array (file_Aj,'i')
    Ap = read_binary_array (file_Ap,'i')

    A = sparse.csr_matrix ((Ax,Aj,Ap))
    
    plt.spy(A, precision=0, marker='.', aspect='equal')
    plt.show()

# end of plot_csr_matrix

if __name__ == "__main__":

    if len(sys.argv) != 4:
        usage()
        sys.exit(1)

    # Read from binary file the CSR Matrix
    fx = sys.argv[1]
    fj = sys.argv[2]
    fp = sys.argv[3]

    print fx, fj, fp

    # Plot sparse CSR matrix
    plot_csr_matrix(fx, fj, fp)





