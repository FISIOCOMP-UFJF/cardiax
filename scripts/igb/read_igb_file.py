#!/usr/bin/env python

"""
Function to read a binary .igb.gz file

Example: [vm, hd] = read_igb_slice(filename)   
Input  : filename
Output : vm array
         hd header

Bernardo M. Rocha
"""

import sys, os, pdb
import numpy, stat

from numpy import *
from scipy.io.npfile import npfile
#from scipy.io.fopen import fopen
from igb_header import igb_header

def read_igb_header(fileName):
      
    try:  
        fd  = open(fileName,'r')
    except IOError:
        print " IOError: File %s not found" % fileName
        exit(-1)
        
    hdl = fd.read(1024)
    aux = hdl.split()
    str = [each.split(':') for each in aux]
    
    # create header
    hd = igb_header()
    hd.initialize()
    
    for tok in str:
        key = tok[0]
        val = tok[1]

        if   key == 'x':        hd.x        = int(val)
        elif key == 'y':        hd.y        = int(val)
        elif key == 'z':        hd.z        = int(val)
        elif key == 't':        hd.t        = int(val)
        elif key == 'dim_x':    hd.dim_x    = float(val)
        elif key == 'dim_y':    hd.dim_y    = float(val)
        elif key == 'dim_z':    hd.dim_z    = float(val)
        elif key == 'dim_z':    hd.dim_z    = float(val)
        elif key == 'dim_t':    hd.dim_t    = float(val)
        elif key == 'unites_x': hd.unites_x = val
        elif key == 'unites_y': hd.unites_y = val
        elif key == 'unites_z': hd.unites_z = val
        elif key == 'unites_t': hd.unites_t = val
        elif key == 'unites':   hd.unites   = val
        elif key == 'type':     hd.type     = val
        elif key == 'systeme':  hd.systeme  = val
        elif key == 'fac_t':    hd.fac_t    = float(val)
        elif key == 'facteur':  hd.facteur  = float(val)
        elif key == 'zero':     hd.zero     = float(val)
        else:
            print " ERROR: Unrecognized token " + key + " in read_ig_header()\n"
            exit(-1)
            
    fd.close()
    return hd

def read_igb_slice(filename):
     
    # get the header from the binary file
    hd = read_igb_header(filename)
    
    # obtain the size of filename
    filestats = os.stat(filename)
    filesize  = filestats[stat.ST_SIZE]

    if hd.systeme == 'big_endian':
        fopenstr='ieee-be'
    else:
        fopenstr='ieee-le'

    # setup time slices to be read
    t_slices = xrange(1,hd.t+1)
    if size(t_slices) == 0:
        t_slices = 1
        print " WARNING: Trying to read at least one time slice!"
    
    # how many time slices we are going to read?
    n_slices = len(t_slices)
    
    # expected data size
    data_in_file = hd.x * hd.y * hd.z * hd.t
    
    # data type
    if   hd.type == 'float':  dbytes = 4; dtype = 'float32';
    elif hd.type == 'double': dbytes = 8; dtype = 'float64';
    else:                     dbytes = 4; dtype = '';

    # check if file is big enough
    if filesize - 1024 < data_in_file * dbytes:
        # adjust data_size
        dtoks = ( filesize - 1024 ) / dbytes
        possible_t_slices = floor( dtoks / hd.x / hd.y / hd.z )
        data_in_file = hd.x * hd.y * hd.z * possible_t_slices

        msg_l1 = 'Mismatch file size with file header: \n'
        msg_l2 = 'Less time slices will be read (%d instead of %d)' % (possible_t_slices, hd.t)
        print " WARNING: " + msg_l1 + msg_l2

    
    # ieee-be' or 'b' IEEE floating point with big-endian byte ordering
    #'ieee-le' or 'l' IEEE floating point with little-endian byte ordering
    if   fopenstr == 'ieee-le':
        fmt = 'l'
    elif fopenstr == 'ieee-be':
        fmt = 'B'
    else:
        print " ERROR in read_igb_slice: specified format doesn't exist!"
    
    # open file in big or little endian machine format    
    #fh = fopen (filename,'r',fmt)
    fh = npfile (filename,'r',fmt)
    
    # skip 1024 header bytes
    #header = fh.read(1024,'char')
    dt = numpy.dtype('uint8')
    header = fh.read_array(dt, 1024)
       
    # preallocate memory
    if    hd.type == 'float':  data = zeros( (hd.x, hd.y, hd.z, n_slices) )
    else:                      data = zeros( (hd.x, hd.y, hd.z, n_slices) )
    
    # size of one time slice
    slice_size = hd.x * hd.y * hd.z
    
    # read data till end of file
    actual_timesteps = 0
    for i in xrange(n_slices):
        # compute position of time slice i
        pos = (t_slices[i] - 1) * slice_size * dbytes + 1024
        fh.seek(pos)
        
        #slice_buf = fh.read(slice_size, dtype)
        
        dt = numpy.dtype(dtype)
        slice_buf = fh.read_array(dt, slice_size)
        count     = size(slice_buf)
        
        if count == slice_size:
            #print " read_igb_slices: Reading time step %d of %d " % (t_slices[i],hd.t)
            data[:,:,:,i] = slice_buf.reshape(hd.x, hd.y, hd.z)
            actual_timesteps = actual_timesteps + 1
        else:
            print " read_igb_slices: Incomplete time step %d of %d " % (i,n_slices)
    
    hd.t = actual_timesteps
    fh.close()
    
    return data, hd