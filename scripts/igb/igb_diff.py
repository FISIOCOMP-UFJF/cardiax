#!/usr/bin/env python

import sys, gzip, pdb
from numpy import zeros
from igb_header import igb_header
from igb_read import *
from igb_write import *

def igb_diff(igbfile1, igbfile2):
    
    [vm1, hd1] = read_igb_slice(igbfile1)
    [vm2, hd2] = read_igb_slice(igbfile2)
    
    Vm1 = vm1.squeeze()
    Vm2 = vm2.squeeze()
    
    print "%s header" % igbfile1
    print hd1
    print "%s header" % igbfile2
    print hd2
    
    shp = shape(Vm1)
    Vm_diff = zeros(shp)
            
    for i in xrange(shp[1]):
        Vm_diff[:,i] = Vm1[:,i] - Vm2[:,i]
    
    write_igb (Vm_diff, hd1, 'vm_diff.igb', gzipped=False)
    print "Vm_diff was written successfully!"
    
    pdb.set_trace()
    
    return Vm_diff

if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print "Usage: igb_diff [igb_file_1] [igb_file_2]"
        exit(-1)
    
    igb_diff(sys.argv[1], sys.argv[2])
                
    #pdb.set_trace()


