#!/usr/bin/env python
#
# Functions for computing the fiber, sheet and sheet-normal
# directions for an ellipsoid mesh of the left ventricle.
# Bernardo M. Rocha
#

import sys, time, os
import numpy as np
from math import *
from scipy.optimize import fsolve
from fiberUtil import rotateAroundAxis

def calcfib(pts):
    """
    Given an ellipsoid mesh file in Medit format (.mesh), read the points
    and then computes the (r,theta,phi) coordinates, compute the local
    basis vectors e_r, e_theta and e_phi and finally computes the fiber 
    orientation a(x) at each node.

    r=0 endo
    r=1 epi
    
    Ref: Spread of Excitation in 3D Models of Anisotropic Cardiac Tissue
         II. Effects of Fibers Architecture and Ventricular Geometry
         P. C. Franzone et al.
    """
    num = np.shape(pts)[0]
    vx,vy,vz = np.zeros((num)),np.zeros((num)),np.zeros((num))
    rvec,tvec,pvec = np.zeros((num)),np.zeros((num)),np.zeros((num))

    for i in range(num):
        x,y,z = pts[i,0], pts[i,1], pts[i,2]
        vx[i], vy[i], vz[i] = x, y, z

        #a1, a2 = 1.1, 1.7
        #c1, c2 = 2.5, 2.8
        
        #a1, a2 = 28000, 43000
        #c1, c2 = 64000, 70000

        #a1, a2 = 19000, 29000
        #c1, c2 = 42000, 47000            

        # 2x menor
        #a1, a2 = 0.5*1.8e+4, 0.5*2.9e+4
        #c1, c2 = 0.5*4.2e+4, 0.5*4.7e+4

        # original
        a1, a2 = 1.8e+4, 2.9e+4
        c1, c2 = 4.2e+4, 4.7e+4       

        a = lambda r: a1 + r*(a2-a1)
        c = lambda r: c1 + r*(c2-c1)
        F = lambda u: [a(u[0])*cos(u[1])*cos(u[2]) - x, 
                       a(u[0])*cos(u[1])*sin(u[2]) - y, 
                       c(u[0])*sin(u[1])           - z]

        # initial guess
        r0 = 0.0                  # sqrt(x*x+y*y+z*z)
        t0 = 0.0                  # acos(z/r)
        p0 = atan2(y,x)           # ok
        v0 = np.array([r0,t0,p0])

        # solve nonlinear problem
        sol = fsolve(F, v0, xtol=1.0e-8)
        
        # save
        rvec[i],tvec[i],pvec[i] = sol[0],sol[1], sol[2]

    print("")
    print("information")
    print("===========")    
    print("r     [%f, %f] " % (rvec.min(),rvec.max()))
    print("theta [%f, %f] " % (tvec.min(),tvec.max()))
    print("phi   [%f, %f] " % (pvec.min(),pvec.max()))
    print("x     [%f, %f] " % (vx.min(),vx.max()))
    print("y     [%f, %f] " % (vy.min(),vy.max()))
    print("z     [%f, %f] " % (vz.min(),vz.max()))

    fib = open("_output/fibras_er.txt","w")
    for i in range(num):
        r,t,p = rvec[i],tvec[i],pvec[i]
        er = np.array([cos(t)*cos(p), cos(t)*sin(p), sin(t)])
        fib.write("%d %f %f %f 0 0 0 0 0 0\n" % (i,er[0],er[1],er[2]))
    fib.close()

    fib = open("_output/fibras_et.txt","w")
    for i in range(num):
        r,t,p = rvec[i],tvec[i],pvec[i]
        et = np.array([-r*sin(t)*cos(p), -r*sin(t)*sin(p), r*cos(t)])
        fib.write("%d %f %f %f 0 0 0 0 0 0\n" % (i,et[0],et[1],et[2]))
    fib.close()

    fib = open("_output/fibras_ep.txt","w")
    for i in range(num):
        r,t,p = rvec[i],tvec[i],pvec[i]
        ep = np.array([-r*cos(t)*sin(p), r*cos(t)*cos(p), 0.0])
        fib.write("%d %f %f %f 0 0 0 0 0 0\n" % (i,ep[0],ep[1],ep[2]))
    fib.close()

    fer = open("_output/fibras_base_er.txt","w")
    fet = open("_output/fibras_base_et.txt","w")
    fep = open("_output/fibras_base_ep.txt","w")
    fib = open("_output/fibras_a.txt","w")

    fib.write('%d\n' % num)

    for i in range(num):
        r,t,p = rvec[i],tvec[i],pvec[i]

        #print(r)
        if (r<0.0):
            r = 0.000001
        elif (r>1.0):
            r = 1.0

        # local base vectors
        er = np.array([cos(t)*cos(p), cos(t)*sin(p), sin(t)])
        et = np.array([-r*sin(t)*cos(p), -r*sin(t)*sin(p), -r*cos(t)])
        ep = np.array([-r*cos(t)*sin(p), r*cos(t)*cos(p), 0.0])

        # normalize        
        ver = er/np.linalg.norm(er)
        vet = et/np.linalg.norm(et)
        vep = ep/np.linalg.norm(ep)

        # compute fiber direction
        #alfa = lambda rr: ((2.0*pi)/3.0)*(1.0-rr) - (pi/4.0)
        #alfa = lambda rr: (2.0*pi*(1.0-rr))/3.0 - (pi/4.0)
        alfa = lambda rr: (-(2.0*pi)/3.0)*(1.0-rr) + (pi/3.0)
        #alfa = lambda x: ((pi*1.0)/1.0)*(1.0 - x) - (pi/4.0)
        #alfa = lambda rr: ((pi*205.0)/180.0)*(1.0 - rr) - (70*pi/180.0)

        a = vep * cos(alfa(r)) + vet * sin(alfa(r))
        a = a/np.linalg.norm(a)
        
        # write vectors to files
        fer.write("%d %f %f %f 0 0 0 0 0 0\n" % (i+1,ver[0],ver[1],ver[2]))
        fet.write("%d %f %f %f 0 0 0 0 0 0\n" % (i+1,vet[0],vet[1],vet[2]))
        fep.write("%d %f %f %f 0 0 0 0 0 0\n" % (i+1,vep[0],vep[1],vep[2]))
        fib.write("%f %f %f \n"  % (a[0],  a[1],  a[2]))
        #fib.write("%f %f %f "  % (-1, -1, -1))
        #fib.write("%f %f %f\n" % (-1, -1, -1))

    fib.close()
    fer.close()
    fet.close()
    fep.close()     

    print("\ndone")
    # end of calcfib

def calcSheet(pts):
    """
    We use s_0 = ( x/sqrt(x^2 + y^2) , y/sqrt(x^2 + y^2), 0) as an 
    approximation for the sheets direction. The fiber direction f0 is
    defined orthogonal to s0 in the planes z = constant; and it is 
    rotated by an angle of 45 deg with respect to the sheet axis 
    to get an approximate direction on the outer surface.
    
    Ref: Orthotropic active strain models for the numerical simulation
         of cardiac biomechanics, Simone Rossi et al, IJBNME, 2013.
    """
    num = np.shape(pts)[0]
    vx,  vy,  vz   = np.zeros((num)),np.zeros((num)),np.zeros((num))
    rvec,tvec,pvec = np.zeros((num)),np.zeros((num)),np.zeros((num))

    for i in range(num):
        x,y,z = pts[i,0], pts[i,1], pts[i,2]
        vx[i], vy[i], vz[i] = x, y, z

        #a1, a2 = 1.1, 1.7
        #c1, c2 = 2.5, 2.8

        #a1, a2 = 28000, 43000
        #c1, c2 = 64000, 70000

        # 10x menor
        #a1, a2 = 1800, 2900
        #c1, c2 = 4200, 4700

        # 2x menor
        #a1, a2 = 0.5*1.8e+4, 0.5*2.9e+4
        #c1, c2 = 0.5*4.2e+4, 0.5*4.7e+4

        a1, a2 = 0.5 * 1.8e+4, 0.5 * 2.9e+4
        c1, c2 = 0.5 * 4.2e+4, 0.5 * 4.7e+4       

        # rabbit heart
        #a1, a2 =  6900, 12700
        #c1, c2 = 18000, 19300

        a = lambda r: a1 + r*(a2-a1)
        c = lambda r: c1 + r*(c2-c1)

        F = lambda u: [a(u[0])*cos(u[1])*cos(u[2]) - x, 
                       a(u[0])*cos(u[1])*sin(u[2]) - y, 
                       c(u[0])*sin(u[1])           - z]

        # initial guess
        r0 = 0.0                  # sqrt(x*x+y*y+z*z)
        t0 = 0.0                  # acos(z/r)
        p0 = atan2(y,x)           # ok
        v0 = np.array([r0,t0,p0])

        # solve nonlinear problem
        sol = fsolve(F, v0, xtol=1.0e-9)

        rvec[i],tvec[i],pvec[i] = sol[0],sol[1],sol[2]

    print("r     [%f, %f] " % (rvec.min(),rvec.max()))
    print("theta [%f, %f] " % (tvec.min(),tvec.max()))
    print("phi   [%f, %f] " % (pvec.min(),pvec.max()))
    print("x     [%f, %f] " % (vx.min(),vx.max()))
    print("y     [%f, %f] " % (vy.min(),vy.max()))
    print("z     [%f, %f] " % (vz.min(),vz.max()))

    # f0, s0 and n0 directions
    ff = open("output_fiber_f0.txt","w")
    fs = open("output_sheet_s0.txt","w")
    fn = open("output_snorm_n0.txt","w")
    fo = open("output_fsn.txt","w")
    fo.write("0\n")
    ff.write('%d\n' % num)
    fs.write('%d\n' % num)
    fn.write('%d\n' % num)

    for i in range(num):
        x,y,z = vx[i], vy[i], vz[i]

        # sheet direction
        rxy = sqrt(x*x + y*y)
        s = np.array([x/rxy, y/rxy, 0.0])

        # fiber direction
        f = np.array([s[1], -s[0], s[2]]) #-s[0]
        r = rvec[i] #abs(rvec[i])

        # angle to rotate f
        #alfa = lambda rr: ((pi*140)/180.0)*(1.0 - rr) - (70*pi/180.0)
        #alfa = lambda rr: ((pi*90.0)/180.0)*(1.0 - rr) - (45*pi/180.0)
        # 70
        alfa = lambda rr: -((pi*145.0)/180.0)*(1.0 - rr) + (85*pi/180.0)

        f = rotateAroundAxis(f, s, alfa(r))

        # normal - cross product 
        n = np.cross(f,s)

        # check
        #print("orthogonals: fs, fn, sn ?", np.dot(f,s), np.dot(f,n), np.dot(s,n))

        f = f/np.linalg.norm(f)
        s = s/np.linalg.norm(s)
        n = n/np.linalg.norm(n)

        ff.write("%f %f %f\n" % (f[0],f[1],f[2]))
        fs.write("%f %f %f\n" % (s[0],s[1],s[2]))      
        fn.write("%f %f %f\n" % (n[0],n[1],n[2]))
        fo.write("%f %f %f "  % (f[0],f[1],f[2]))
        fo.write("%f %f %f "  % (s[0],s[1],s[2]))
        fo.write("%f %f %f\n" % (n[0],n[1],n[2]))

    ff.close()
    fs.close()
    fn.close()
    fo.close()
    print("Done")

if __name__ == "__main__":

    if (len(sys.argv) < 2):
        print("\n Usage: ellipsoidFibGen <file.mesh> OR ellipsoidFibGen <basenameCarp>\n")
        sys.exit(-1)

    # parse and check input
    infile = sys.argv[1]

    # extract nodes of the mesh
    infile = infile
    if(infile[-4:]=='.msh'):
        infile = infile
    else:
        infile = infile + '.pts'

    if (not os.path.isfile(infile)):
        print("\n Error: the input %s does not exist.\n" % (infile))
        sys.exit(-1)

    print('Reading file %s' % infile)
    finp = open(infile)
    npts = int(finp.readline())
    vpts = np.zeros((npts,3))
    for i in range(npts):
        line = finp.readline()
        x,y,z = map(float,line.split())
        vpts[i,0], vpts[i,1], vpts[i,2] = x, y, z    
    finp.close()

    # compute fibers
    #print("Computing fiber directions - first fiber, sheet and then normal")
    #calcfib(vpts)

    # compute sheet
    print("Computing fiber directions - first sheet, fiber and then normal")
    calcSheet(vpts)
    print("Vectors defined at the nodes!")

# end of main
