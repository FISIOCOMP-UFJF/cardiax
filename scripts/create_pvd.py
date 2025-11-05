#!/usr/bin/env python
# coding: utf-8
#
# Bernardo M. Rocha, Fev, 2014
#

import sys
from os import listdir
from os.path import isfile, join

def file_num_key(item):
	last = item.split('_')[-1]
	numb = last.split('.')[0]
	return int(numb)

if __name__ == "__main__":

	if (len(sys.argv) != 3):
		print("\n Usage: create_pvd.py [directory with vtu files] [outputname.pvd]\n")
		sys.exit(1)

	mypath  = sys.argv[1]
	pvdname = sys.argv[2]

	myfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
	myfiles = sorted(myfiles, key=file_num_key)

	fp = open(pvdname,"w")
	fp.write("<?xml version=\"1.0\"?>\n")
	fp.write("<VTKFile type=\"Collection\" version=\"0.1\" ")
	fp.write("byte_order=\"LittleEndian\" ")
	fp.write("compressor=\"vtkZLibDataCompressor\">\n")
	fp.write("  <Collection>\n")

	for i in range(len(myfiles)):
		outname = mypath + myfiles[i]
		fp.write("<DataSet timestep=\"")
		fp.write("%d \" group=\"\" part=\"0\" " % i)
		fp.write("file=\"%s\"/>\n" % outname)

	fp.write("  </Collection>\n")
	fp.write("</VTKFile>\n")
	fp.close()
