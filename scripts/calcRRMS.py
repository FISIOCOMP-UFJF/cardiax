#!/usr/bin/env python

"""
Python module to compute the RRMS (relative root-mean-square) norm error
 - calc_rrms: for identical meshes
 - calc_map_rrms: for different meshes (compute at the common nodes only)
 
Bernardo M. Rocha, 2009
"""

import sys, pdb
import numpy as np
from igb.igb_read import *
from sctools import read_map_array

def calc_rrms (ve_name, vh_name):

  print "Computing RRMS"

  ve, vh = load_data(ve_name, vh_name)

  # reshape arrays to 1d arrays
  veshp = np.shape(ve)
  #ve = ve.reshape((veshp[0], veshp[1]*veshp[2]))

  vhshp = np.shape(vh)
  #vh = vh.reshape((vhshp[0], vhshp[1]*vhshp[2]))

  shp = np.shape(vh)
  time_steps = shp[0]

  # compute RRMS norm - call external module
  rrms(time_steps, ve, vh)

# end of calc_rrms

def rrms(time_steps, ve, vh):
  """
  This function is not optmized/implemented in Cython, uses normal Python/numpy
  """
  aux1 = 0.0
  aux2 = 0.0
  t1 = 0.0
  t2 = 0.0

  for t in xrange(time_steps):
    # slice in time t
    se = ve[t,:]
    sh = vh[t,:]
    aux1 += np.sum( (se-sh)*(se-sh) )
    aux2 += np.sum( se*se )
    t1 = np.sum( (se-sh)*(se-sh) )
    t2 = np.sum( se*se )

    err1 = 100.0 * ( np.sqrt(t1) / np.sqrt(t2)	)
    print "RRMS error at %d ms : %f" % (t,err1)

  err = 100.0 * ( np.sqrt(aux1) / np.sqrt(aux2)	)
  print "RRMS error over time: ", err, "%"

# end of rrms

def load_data(ve_name, vh_name):

  # open exact Vm IGB
  if os.path.isfile(ve_name):
    [ve, veh] = read_igb_slice(ve_name, is_gzipped=True)
    ve = ve.squeeze()
  else:
    print "Error: the file %s does not exist." % ve_name
    sys.exit(1)

  # open aprox Vm IGB
  if os.path.isfile(vh_name):
    [vh, vhh] = read_igb_slice(vh_name, is_gzipped=True)
    vh = vh.squeeze()
  else:
    print "Error: the file %s does not exist." % vh_name
    sys.exit(1)

  # show information
  print np.shape(ve)
  print np.shape(vh)

  return ve, vh

# end of load_data

if __name__ == "__main__":

  if len(sys.argv) < 3:
    print "\n Usage: calcRRMS.py <vm_hyb> <vm_tet> <map_hyb_tet>\n"; exit(1);

  ven = sys.argv[1] # nome da solucao exata
  vhn = sys.argv[2] # nome da solucao aprox

  # compute RRMS with map file
  if len(sys.argv) > 3:
    map_name = sys.argv[3]
    calc_map_rrms(ven, vhn, map_name)
  else:
    calc_rrms (ven, vhn)

# end of main
