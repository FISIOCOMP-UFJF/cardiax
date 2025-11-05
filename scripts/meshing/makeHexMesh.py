#!/usr/bin/env python

import sys, pdb
import numpy as np

def make_hex_mesh(basename, num_x, num_y, num_z):

    ptsname = basename + ".pts"
    bdnname = basename + ".bdn"
    elename = basename + ".elem"

    fp = open(ptsname, 'w')
    fb = open(bdnname, 'w')
    fe = open(elename, 'w')

    nx = int(num_x)
    ny = int(num_y)
    nz = int(num_z)

    sx, sy, sz = 10000.0, 100.0, 100.0

    dx = sx/nx
    dy = sy/ny
    dz = sz/nz

    x = np.arange(0.0,sx+dx,dx)
    y = np.arange(0.0,sy+dy,dy)
    z = np.arange(0.0,sz+dz,dz)

    numNodes = len(x) * len(y) * len(z)
    XYZ = np.zeros((numNodes,3))

    fp.write("%d\n" % numNodes)

    c = 0
    for k in range(len(z)):
        for i in range(len(x)):
            for j in range(len(y)):
                XYZ[c][0] = x[i]
                XYZ[c][1] = y[j]
                XYZ[c][2] = z[k]
                # trick for flagshyp
                #code = 0
                #if (x[k]==0.0 and y[j]==0.0): code = 3
                #elif (x[k]==0.0): code = 1
                #elif (y[j]==0.0): code = 2
                #fp.write("%d %d %f %f %f\n" % (c+1,code,x[k],y[j],z[i]))
                fp.write("%f %f %f\n" % (x[i],y[j],z[k]))
                c = c + 1

    # build connectivity
    numElems = nx*ny*nz
    fe.write("%d\n" % numElems)

    c = 0
    Nx,Ny,Nz = nx+1,ny+1,nz+1
    for k in range(nz):
        for i in range(nx): 
            for j in range(ny):
                # trick for flagshyp
                #fe.write("%d 1 %d %d %d %d %d %d %d %d\n" % (c+1,
                #    i*Nx+j+k*(Ny*Nz)+1, i*Nx+j+1+k*(Ny*Nz)+1, 
                #    (i+1)*Nx+j+1+k*(Ny*Nz)+1, (i+1)*Nx+j+k*(Ny*Nz)+1,
                #    i*Nx+j+(k+1)*(Ny*Nz)+1, i*Nx+j+1+(k+1)*(Ny*Nz)+1, 
                #    (i+1)*Nx+j+1+(k+1)*(Ny*Nz)+1, (i+1)*Nx+j+(k+1)*(Ny*Nz)+1))
                fe.write("Hx %d %d %d %d %d %d %d %d 1\n" % (
                    i*Nx+j+k*(Ny*Nz),
                    i*Nx+j+1+k*(Ny*Nz), 
                    (i+1)*Nx+j+k*(Ny*Nz),
                    (i+1)*Nx+j+1+k*(Ny*Nz),
                    i*Nx+j+(k+1)*(Ny*Nz),
                    i*Nx+j+1+(k+1)*(Ny*Nz), 
                    (i+1)*Nx+j+(k+1)*(Ny*Nz),
                    (i+1)*Nx+j+1+(k+1)*(Ny*Nz)))
                c = c + 1
    
    # boundary conditions
    c = 0
    bdlst=[]
    for i in range(numNodes):
        px,py,pz = XYZ[i][0],XYZ[i][1],XYZ[i][2]
        if (px == 0.0 or px == 1.0 or 
            py == 0.0 or py == 1.0 or 
            pz == 0.0 or pz == 1.0):
            bdlst.append( (i,0.0) )
            c = c + 1
    numConsNodes = c
    fb.write("%d\n" % numConsNodes)
    for i in range(numConsNodes):
        fb.write("%d\t%f\n" % (bdlst[i][0], bdlst[i][1]))

    # fibers
    
    
    fp.close()
    fb.close()
    fe.close()
    print "Done."
    
 
if __name__ == "__main__":

    basename = sys.argv[1]
    in_args = sys.argv[2:]
    in_args_len = len(in_args)
    
    make_hex_mesh (basename, in_args[0], in_args[1], in_args[2])

