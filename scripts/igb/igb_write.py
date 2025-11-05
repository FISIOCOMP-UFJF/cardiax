#!/usr/bin/env python

import sys, os, time, gzip, stat, string, pdb
from numpy import zeros
from igb_header import igb_header
from read_igb_file import *

"""
IGB write file
Bernardo M. Rocha
"""

def write_igb(data, hd, filename, gzipped=False):
    
    t0 = time.time()
    
    if data.ndim == 2:                   # if data is (nodes, time_slices)
        data = data[:,newaxis,newaxis,:]
    
    if data.ndim == 3:
        data = data[:,:,:,newaxis]

    n_dims = data.ndim
    d_data = zeros((data.ndim))
    for i in xrange(data.ndim):
        d_data[i] = data.shape[i]

    t1 = time.time()
    pData = pack_float_data(data.squeeze())
    t2 = time.time()

    print " packFloatData: %.5f seconds wall time" % (t2-t1)

    if d_data[0] == hd.x and d_data[1] == hd.y and d_data[2] == hd.z and d_data[3] == hd.t:
        
        if gzipped:
            gzfilename = '%s.gz' % filename
            fh = gzip.GzipFile(filename=gzfilename,mode='wb')
        else:
            fh = open (filename,'wb')

        write_header(hd, fh)
        
        fh.write(pData)
        fh.close()
    else:
        print " ERROR: Header-Data mismatch\n"

    print " write_igb    : %.5f seconds wall time" % (time.time() - t0)
    
def write_header(hd, fh):
    """
    write the 1024 bytes IGB format header (hd) in file (fh)
    """    
    # header size in bytes
    hs = 1024

    # create write buffer
    buffer = [ ' ' for i in xrange(hs-1)]

    # create header substrings
    N  = 'x:%d y:%d z:%d t:%d' % (hd.x, hd.y,hd.z,hd.t)
    DT = 'type:%s' % (hd.type)
    SY = 'systeme:%s \r\n' % (hd.systeme)       
    D  = ''
    if (not isnan(hd.dim_x)): D += 'dim_x:%f ' % hd.dim_x
    if (not isnan(hd.dim_y)): D += 'dim_y:%f ' % hd.dim_y
    if (not isnan(hd.dim_z)): D += 'dim_z:%f ' % hd.dim_z
    if (not isnan(hd.dim_t)): D += 'dim_t:%.2f \r\n' % hd.dim_t
    #D  = 'dim_x:%f dim_y:%f dim_z:%f dim_t:%f \r\n' % (hd.dim_x, hd.dim_y, hd.dim_z, hd.dim_t)
    
    U = ''
    if (hd.unites_x != ''): U += 'unites_x:%s ' % hd.unites_x
    if (hd.unites_y != ''): U += 'unites_y:%s ' % hd.unites_y
    if (hd.unites_z != ''): U += 'unites_z:%s ' % hd.unites_z
    if (hd.unites_t != ''): U += 'unites_t:%s ' % hd.unites_t
    if (hd.unites   != ''): U += 'unites:%s \r\n'   % hd.unites
            
    #U  = 'unites_x:%s unites_y:%s unites_z:%s unites_t:%s unites:%s \n' % (hd.unites_x, hd.unites_y, hd.unites_z, hd.unites_t, hd.unites)
    F  = 'facteur:%d' % (hd.facteur)
    Z  = 'zero:%d \r\n' % (hd.zero)

    # concatenate header substrings
    hd_str = '%s %s %s %s %s %s %s' % (N, DT, SY, D, U, F, Z)
        
    # fill rest with empty lines, each line not longer than 80 chars
    lines = int(floor( (1023 - len(hd_str) ) / 78 ))
   
    ln = [' ' for i in xrange(75)]    
    ln_str = string.join(ln, '')
    
    for i in xrange(lines):
        hd_str = '%s %s \r\n' % (hd_str, ln_str)

    # convert to character array and write to file
    for i in xrange(len(hd_str)):
        buffer[i] = hd_str[i]
    buffer_str = string.join(buffer, '')
    fh.write(buffer_str)    

    # write header termination character (seems to be 0x0C)
    hd_term = chr(12)
    fh.write(hd_term)

################################################################################

def isnan(x):
    return str(x) == 'nan'

def pack_float_data(data):
    """
    pack binary data to be written in a binary file in 'little_endian' format
    data must has 2 dimensions (ndim=2)
    data(x,t) -> x node
              -> t time slice
    """
    import struct    
    packedData=''
    size = shape(data)
    for j in xrange(size[1]):
        for i in xrange(size[0]):
            packedData+=struct.pack('<f',data[i,j])
    return packedData

def create_igb(hd, filename, gzipped=True):
 
    if gzipped:
        gzfilename = '%s.gz' % filename
        fh = gzip.GzipFile(filename=gzfilename,mode='wb')
    else:
        fh = open (filename,'wb')

    write_header(hd, fh)
    return fh

def close_igb(fh):   
    fh.close()
    
def append_data_to_igb(data, fh):
    """para 2d"""
    import struct    
    packedData=''
    size = shape(data)
    for i in xrange(size[0]):
        for j in xrange(size[1]):
            packedData+=struct.pack('<f',data[i,j])

    fh.write(packedData)

if __name__ == "__main__":
       
    ### test case for the module
    
    # resolution 100um - UCLA_RAB model    
    igb_file = '/data/sim/space_constant/sc_ucla_long_0100um/vm.igb'
    pts_file = '/data/sim/space_constant/sc_ucla_long_0100um/sc_ucla_long_0100um_i.pts'
    [vm, hd] = read_igb_slice(igb_file)
    Vm = vm.squeeze()
    pdb.set_trace()
    write_igb (Vm, hd, 'teste.igb')

