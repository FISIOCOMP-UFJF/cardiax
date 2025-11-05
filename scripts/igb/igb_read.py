#!/usr/bin/env python

import sys, os, gzip, stat, numpy, struct, time
from numpy import *
from igb_header import igb_header
from scipy.io.numpyio import fwrite, fread

"""
Module for reading IGB files

Example: [vm, hd] = read_igb_slice(filename) where vm is [time, x, y, z]
Input  : filename
Output : vm array
         hd header

Bernardo M. Rocha
"""

def read_igb_header(fileName, is_gzipped):

    try:
        if is_gzipped:
            fd = gzip.GzipFile(filename=fileName,mode='rb')
        else:
            fd = open(fileName,'rb')
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

def read_igb_slice (filename, is_gzipped=False):
    """
    Reads IGB slice. If gzipped it uncompress and output a binary file
    that is used to read the data
    """
    if is_gzipped:
        igbFile = gunzipFile(filename)
    else:
        igbFile = filename

    hd = read_igb_header(igbFile, is_gzipped=False)

    filestats = os.stat(igbFile)
    filesize  = filestats[stat.ST_SIZE]

    if   hd.systeme == 'big_endian'   : byteswap=1
    elif hd.systeme == 'little_endian': byteswap=0

    # setup time slices to be read
    t_slices = xrange(1,hd.t+1)
    if size(t_slices) == 0:
        t_slices = 1
        print " WARNING: Trying to read at least one time slice!"

    # how many time slices we are going to read ?
    n_slices = len(t_slices)

    # expected data size
    data_in_file = hd.x * hd.y * hd.z * hd.t

    ## FORCE MY CASE in case of MEMORY ERROR ##
    #n_slices = 1001

    # data type
    if hd.type == 'float':  
        dbytes = 4
        dtype  = 'f'
        data   = zeros( (n_slices, hd.x, hd.y, hd.z), dtype=float32 )
    elif hd.type == 'double':
        dbytes = 8
        dtype = 'd'
        data  = zeros( (n_slices, hd.x, hd.y, hd.z), dtype=float64 )        

    # open data
    fh = open(igbFile,'rb')
   
    # size of one time slice
    slice_size = hd.x * hd.y * hd.z

    # read data till end of file
    actual_timesteps = 0
    for i in xrange(n_slices):
        # compute position of time slice i
        pos = (t_slices[i] - 1) * slice_size * dbytes + 1024        
        fh.seek(pos)

        slice_buf = fread(fh, slice_size, dtype, dtype, byteswap)

        count = size(slice_buf)

        if count == slice_size:
            data[i,:,:,:] = slice_buf.reshape(hd.x, hd.y, hd.z)
            actual_timesteps = actual_timesteps + 1
        #else:
        #    print " read_igb_slices: Incomplete time step %d of %d " % (i,n_slices)

    hd.t = actual_timesteps
    fh.close() 

    # if gzipped, remove temporary uncompressed file
    if is_gzipped:
        os.remove(igbFile)

    return data, hd

def gunzipFile (filename):
    """
    Uncompress file and create a temporary with the contents of compressed file
    """
    pass
    f = gzip.open(filename)
    t = filename[:-3]
    g = open(t,'wb')
    while 1:
        chunk = f.read(1024)
        if not chunk:
            break
        g.write(chunk)
    f.close()
    g.close()
    return t
