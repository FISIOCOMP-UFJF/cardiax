#!/usr/bin/env python

import sys, time, os
import numpy as np

def febioGetBoundary(febMesh, output):

    f = open(febMesh)
    of = file(output, 'w')

    of.write("material_type CompressibleNeoHookean\n")
    of.write("material_coef 2 3.0 2.0\n")
    of.write("problem_type THREE_DIM\n")
    of.write("num_increments 10\n\n")

    line = f.readline()
    while(line.find("<Boundary>") == -1):
        line = f.readline()
       
    of.write("prescribed_displacement\n")

    count = 0
    strpd = ""

    line = f.readline() 
    while(line.find('fix')>0):
        s,c = read_fix(f)
        strpd = strpd + s
        count = count + c
        line = f.readline() 

    print("Done reading <fix> sections")

    print line
    if (line.find('prescribe')>0):
        s,c = read_prescribe(f)    
        strpd = strpd + s
        count = count + c
        line = f.readline() 
        while(line.find('prescribe')>0):
            s,c = read_prescribe(f)
            strpd = strpd + s
            count = count + c
            line = f.readline() 
        
    """
    print("Reading <prescribe> section")
    line = f.readline() # <prescribe>
    line = f.readline()
    while(line.find("</prescribe>") == -1):
        temp = line.split()       
        
        print temp, len(temp)

        node_id = int(temp[1].split('"')[1])
        node_id = node_id - 1
        ftype = temp[2].split('"')[1]             

        if(len(temp)==4):            
            # <node id="296" bc="xyz"/> 1 </node>
            lc = temp[3].split('>')[1]
            lc = lc.split('<')
            pd = float(lc[0])
        else:
            # <node id="296" bc="xyz"/>
            lc = temp[2].split('=')[1]
            pd = 0.0
        
        if (ftype == 'xyz'):
            strpd += '%d 0 %f\n' % (node_id,pd)
            strpd += '%d 1 %f\n' % (node_id,pd)
            strpd += '%d 2 %f\n' % (node_id,pd)
        elif (ftype == 'yz'):
            strpd += '%d 1 %f\n' % (node_id,pd)
            strpd += '%d 2 %f\n' % (node_id,pd)
        elif (ftype == 'xz'):
            strpd += '%d 0 %f\n' % (node_id,pd)
            strpd += '%d 2 %f\n' % (node_id,pd)
        elif (ftype == 'xy'):
            strpd += '%d 0 %f\n' % (node_id,pd)
            strpd += '%d 1 %f\n' % (node_id,pd)
        elif (ftype == 'x'):
            strpd += '%d 0 %f\n' % (node_id,pd)
        elif (ftype == 'y'):
            strpd += '%d 1 %f\n' % (node_id,pd)
        elif (ftype == 'z'):
            strpd += '%d 2 %f\n' % (node_id,pd)
        else:
            print("Error: bc type unknown")
            sys.exit(-1)    

        line = f.readline()
        count = count + 1
    """

    of.write("%d\n" % count)
    of.write("%s" % strpd)
    of.write("end_prescribed_displacement\n")        
    of.close()
    f.close()

def read_prescribe(f):
    print("Reading <prescribe> section")

    count = 0
    strpd = ""

    #line = f.readline() # <prescribe>

    line = f.readline()
    while(line.find("</prescribe>") == -1):
        temp = line.split()       
        node_id = int(temp[1].split('"')[1])
        node_id = node_id - 1
        ftype = temp[2].split('"')[1]             
        
        if(len(temp)==4):            
            # <node id="296" bc="xyz"/> 1 </node>
            lc = temp[3].split('>')[1]
            lc = lc.split('<')
            pd = float(lc[0])            

        if (ftype == 'xyz'):
            strpd += '%d 0 %f\n' % (node_id,pd)
            strpd += '%d 1 %f\n' % (node_id,pd)
            strpd += '%d 2 %f\n' % (node_id,pd)
        elif (ftype == 'yz'):
            strpd += '%d 1 %f\n' % (node_id,pd)
            strpd += '%d 2 %f\n' % (node_id,pd)
        elif (ftype == 'xz'):
            strpd += '%d 0 %f\n' % (node_id,pd)
            strpd += '%d 2 %f\n' % (node_id,pd)
        elif (ftype == 'xy'):
            strpd += '%d 0 %f\n' % (node_id,pd)
            strpd += '%d 1 %f\n' % (node_id,pd)
        elif (ftype == 'x'):
            strpd += '%d 0 %f\n' % (node_id,pd)
        elif (ftype == 'y'):
            strpd += '%d 1 %f\n' % (node_id,pd)
        elif (ftype == 'z'):
            strpd += '%d 2 %f\n' % (node_id,pd)
        else:
            print("Error: bc %s  type unknown" % ftype)
            sys.exit(-1)    

        line = f.readline()
        count = count + 1

    return strpd, count

def read_fix(f):
    print("Reading <fix> section")

    count = 0
    strpd = ""

    line = f.readline()
    print line
    while(line.find("</fix>") == -1):
        temp = line.split()       
        node_id = int(temp[1].split('"')[1])
        node_id = node_id - 1
        ftype = temp[2].split('"')[1]

        if (ftype == 'xyz'):
            strpd += '%d 0 0.0\n' % node_id
            strpd += '%d 1 0.0\n' % node_id
            strpd += '%d 2 0.0\n' % node_id
        elif (ftype == 'yz'):
            strpd += '%d 1 0.0\n' % node_id
            strpd += '%d 2 0.0\n' % node_id
        elif (ftype == 'xz'):
            strpd += '%d 0 0.0\n' % node_id
            strpd += '%d 2 0.0\n' % node_id
        elif (ftype == 'xy'):
            strpd += '%d 0 0.0\n' % node_id
            strpd += '%d 1 0.0\n' % node_id
        elif (ftype == 'x'):
            strpd += '%d 0 0.0\n' % node_id
        elif (ftype == 'y'):
            strpd += '%d 1 0.0\n' % node_id
        elif (ftype == 'z'):
            strpd += '%d 2 0.0\n' % node_id
        else:
            print("Error: bc type unknown")
            sys.exit(-1)                    
            
        line = f.readline()
        count = count + 1

    return strpd, count
# end


if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: febio2gmsh <file.feb> <file.msh>\n")
        sys.exit(-1)

    # parse and check input
    feb_mesh = sys.argv[1]
    out_mesh = sys.argv[2]

    if (not os.path.isfile(feb_mesh)):
        print("\n Error: the input file %s does not exist.\n" % (feb_mesh))
        sys.exit(-1)

    # convert
    febioGetBoundary(feb_mesh, out_mesh)

    print("Done")

# end of main
