#!/usr/bin/env python

import sys, gzip, pdb, time
from numpy import zeros, array, shape, linalg, arange, max
from pyUtils.igb import igb_read as rigb
from pyUtils.igb import igb_read as wigb
from pyUtils.igb import igb_read as higb
from pyUtils.igb import read_igb_file
from pyUtils.sctools import read_array_pts
from pyUtils.petsc.binary_read import petsc_binary_read
from scipy import sparse
from pylab import *

def igb_mapped_norm(igbfile1, igbfile2, mapfile):
    """
    Compute the L2 Norm of two IGB files based on a map file that
    match the node numbering of the two different IGB files.
    """
    
    [Vm1, hd1] = rigb.read_igb_slice(igbfile1)
    [Vm2, hd2] = rigb.read_igb_slice(igbfile2)    
    
    vm1 = Vm1.squeeze(); vm2 = Vm2.squeeze()
    
    shp1  = shape(vm1)
    shp2  = shape(vm2)
    xdim1 = shp1[1]; time1 = shp1[0];
    xdim2 = shp2[1]; time2 = shp2[0]
    
    print ""

    if mapfile is not None:
        # match nodes in hybrid and tetrahedra mesh
        maparray = read_array_pts(mapfile)
        size = len(maparray[:,0])   
        temp = zeros((time1,size))
        
        for t in xrange(time1):
            #count = 0
            for i in xrange(len(maparray[:,0])):
                ind1 = int(maparray[i][0])
                ind2 = int(maparray[i][1])
                temp[t][i] = vm1[t][ind1] - vm2[t][ind2]
                #count = count + 1
    else:        
        # compute the difference
        temp = zeros((time1,xdim1))
        for t in xrange(time1):
            count = 0
            temp[t][:] = vm1[t][:] - vm2[t][:]
    
    erro = zeros(time1)
    
    print '============ P O I N T W I S E   C O M P A R I S O N ==========\n'
    print 'At time t=%d ms' % (20)    
    print ' Maximum norm at time=20ms : ' , max( temp[20,:]  )
    print ' L2 norm at time=20ms      : ' , linalg.norm( temp[20,:]  )
    print '\nAt time t=%d ms' % (100)    
    print ' Maximum norm at time=100ms: ' , max( temp[100,:] )
    print ' L2 norm at time=100ms     : ' , linalg.norm( temp[100,:] )
    print '\nAt time t=%d ms' % (120)    
    print ' Maximum norm at time=120ms: ' , max( temp[120,:] )
    print ' L2 norm at time=120ms     : ' , linalg.norm( temp[120,:] )
    
    n2chk = 450
    print '\nDifference in arrival time (n=450): ' , chk_dif_in_max_atime(vm1[:,n2chk],vm2[:,n2chk])      
        
    for t in xrange(time1):
        #erro[t] = linalg.norm(temp[t,:])
        erro[t] = max(temp[t,:])
    
    #ax = subplot(111)
    #mytime = arange(0,151)
    #ax.plot(mytime, vm1[:,n2chk], mytime, vm2[:,n2chk])
    #ax.legend(('hyb','tet'), 'upper right')
    #show()

    
    ### new stuff ###
    print '\n============= E R R O R   C O M P A R I S O N =============\n'

    # interpolate solution

    # compute error as a column vector
    t = 10
    e = temp[t][:]

    # read and store CARP lumped mass matrix
    massfile1 = '/data/sim/simulacao_1/hyb_75um/output/MatLabDump_Mi.bin'
    M = petsc_binary_read (massfile1,0)
    size = np.size(M)
    data = M
    rows = np.arange(0,size,1)
    cols = np.arange(0,size,1)
    A = sparse.coo_matrix( (data,(rows,cols)) ,(size, size) )
  
    #aux = dot(e,A*e)
    print '\nAt time t=%d ms' % (t)
    print ' Mean-square-root L2 Norm (normal)     : ', compute_L2_error(e,A)
    print ' Mean-square-root L2 Norm (linalg.norm):' , linalg.norm(e)

    print '\n'
    pdb.set_trace()

def compute_L2_error(e, mass):
    aux = dot(e, mass*e)
    return sqrt(aux)

def chk_dif_in_max_atime(vm1_n, vm2_n):
    dvdt1_at_node = diff(vm1_n, n=1)
    dvdt2_at_node = diff(vm2_n, n=1)    
    val = np.max(dvdt1_at_node) - np.max(dvdt2_at_node)
    return val
    
if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print "\n Usage: igb_norm [igb_file_1] [igb_file_2] [mapfile]\n"
        exit(-1)
    
    mapfile = None
    if len(sys.argv) == 4:
        mapfile = sys.argv[3]
    
    igb_mapped_norm(sys.argv[1], sys.argv[2], mapfile)
        



