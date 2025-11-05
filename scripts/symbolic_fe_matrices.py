#!/usr/bin/env python

from sympy import *

def Emn(dphi,m,n):
    val = integrate( dphi[0][m]*dphi[0][n] + dphi[1][m]*dphi[1][n] , 
                     (x,0,1), (y,0,1))
    return val

def Emn2(dphi,m,n):
    val = integrate( dphi[0][m]*dphi[0][n] + dphi[1][m]*dphi[1][n] , 
                     (x,-1.0,1.0), (y,-1.0,1.0))
    return val

def Emn3D(dphi,m,n):
    val = integrate( dphi[0][m]*dphi[0][n] + 
                     dphi[1][m]*dphi[1][n] +
                     dphi[2][m]*dphi[2][n], 
                     (x,0,1), (y,0,1), (z,0,1))
    return val

def Emn3D2(dphi,m,n):
    val = integrate( dphi[0][m]*dphi[0][n] + 
                     dphi[1][m]*dphi[1][n] +
                     dphi[2][m]*dphi[2][n], 
                     (x,-1.0,1.0), (y,-1.0,1.0), (z,-1.0,1.0))
    return val


if __name__ == "__main__":

    x = Symbol('x')
    y = Symbol('y')

    phi=[]
    phi.append( (1-x)*(1-y) )
    phi.append(   x  *(1-y) )
    phi.append( (1-x)*  y   )
    phi.append(   x  *  y   )

    dphi=[]
    dphix=[]
    dphiy=[]

    for i in range(4):
        dphix.append( diff(phi[i],x) )
        dphiy.append( diff(phi[i],y) )
    
    dphi.append(dphix)
    dphi.append(dphiy)

    print "\nBilinear element stiffness matrix"
    for i in range(4):
        for j in range(4):
            print Emn(dphi,i,j), 
        print ""
    
    # nova
        
    phi=[]
    phi.append( 0.25*(1-x)*(1-y) )
    phi.append( 0.25*(1+x)*(1-y) )
    phi.append( 0.25*(1+x)*(1+y) )
    phi.append( 0.25*(1-x)*(1+y) )

    dphi=[]
    dphix=[]
    dphiy=[]

    for i in range(4):
        dphix.append( diff(phi[i],x) )
        dphiy.append( diff(phi[i],y) )
    
    dphi.append(dphix)
    dphi.append(dphiy)

    print "\nBilinear element stiffness matrix"
    for i in range(4):
        for j in range(4):
            print "%.4f " % Emn2(dphi,i,j), 
        print ""

    # 3D de 0 a 1
    z = Symbol('z')

    phi=[]
    phi.append( (1-x)*(1-y)*(1-z) )
    phi.append( x*(1-y)*(1-z) )
    phi.append( (1-x)*y*(1-z) )
    phi.append( x*y*(1-z) )
    phi.append( (1-x)*(1-y)*z )
    phi.append( x*(1-y)*z )
    phi.append( (1-x)*y*z )
    phi.append( x*y*z )

    dphi=[]
    dphix=[]
    dphiy=[]
    dphiz=[]

    for i in range(8):
        dphix.append( diff(phi[i],x) )
        dphiy.append( diff(phi[i],y) )
        dphiz.append( diff(phi[i],z) )
    
    dphi.append(dphix)
    dphi.append(dphiy)
    dphi.append(dphiz)

    print "\nTrilinear element stiffness matrix"
    for i in range(8):
        for j in range(8):
            print "%.4f" % Emn3D(dphi,i,j), 
        print ""
                     
    # 3D de -1 a 1
    z = Symbol('z')

    phi=[]
    phi.append( 0.125*(1-x)*(1-y)*(1-z) )
    phi.append( 0.125*(1+x)*(1-y)*(1-z) )
    phi.append( 0.125*(1+x)*(1+y)*(1-z) )
    phi.append( 0.125*(1-x)*(1+y)*(1-z) )
    phi.append( 0.125*(1-x)*(1-y)*(1+z) )
    phi.append( 0.125*(1+x)*(1-y)*(1+z) )
    phi.append( 0.125*(1+x)*(1+y)*(1+z) )
    phi.append( 0.125*(1-x)*(1+y)*(1+z) )

    dphi=[]
    dphix=[]
    dphiy=[]
    dphiz=[]

    for i in range(8):
        dphix.append( diff(phi[i],x) )
        dphiy.append( diff(phi[i],y) )
        dphiz.append( diff(phi[i],z) )
    
    dphi.append(dphix)
    dphi.append(dphiy)
    dphi.append(dphiz)

    print "\nTrilinear element stiffness matrix"
    for i in range(8):
        for j in range(8):
            print "%.4f" % Emn3D2(dphi,i,j), 
        print ""
        
    
