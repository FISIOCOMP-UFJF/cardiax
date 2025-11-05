import sys, inspect
import numpy as np
from sympy import *
from hyperelastic_base import C

def rows(A):
    return [A[i,:] for i in range(A.rows)]

def grad(U,X):
    return u.jacobian(X)

def symgrad(U, X):
    return sym(grad(U,X))

def sym(A):
    return 0.5*(A + A.T)

def divergence(F, X):
    assert(F.cols == X.rows)
    return sum(diff(F[i],X[i]) for i in range(X.rows))

def trace(T):
    return T[0,0]+T[1,1]+T[2,2]

def contract(A,B):
    s = 0
    for i in range(3):
        for j in range(3):
            s += A[i,j] * B[i,j]
    return s

# Invariants

def I1(C):
    return trace(C)

def I2(C):
    return 0.5 * ((I1(C)*I1(C)) - contract(C,C))

def I3(C):
    return C.det()

def I1_deviatoric(C):
    C_hat = pow(I3(C),-1.0/3.0) * C
    return trace(C_hat)

def I2_deviatoric(C):
    C_hat = pow(I3(C),-1.0/3.0) * C
    return 0.5 * ((I1(C_hat)*I1(C_hat)) - contract(C_hat,C_hat))

def I3_deviatoric(C):
    C_hat = I3(C)**(-1.0/3.0) * C
    return C_hat.det()

# Voigt notation

def stress_deriv_voigt2D(C):
    """
    Given stress tensor derivatives as a fourth order tensor C, 
    this function returns the 2D representation using Voigt notation.

    indexes = ((0,0,0,0), (0,0,1,1), (0,0,0,1),
               (1,1,0,0), (1,1,1,1), (1,1,0,1),
               (0,1,0,0), (0,1,1,1), (0,1,0,1))
    """
    M = Matrix([ [0,0,0], [0,0,0], [0,0,0] ])    

    M[0,0], M[0,1], M[0,2] = C[0,0,0,0], C[0,0,1,1], C[0,0,0,1]
    M[1,0], M[1,1], M[1,2] = C[1,1,0,0], C[1,1,1,1], C[1,1,0,1]
    M[2,0], M[2,1], M[2,2] = C[0,1,0,0], C[0,1,1,1], C[0,1,0,1]

    return M

def stress_deriv_voigt3D(C):
    """
    Given stress tensor derivatives as a fourth order tensor C, 
    this function returns the 3D representation using Voigt notation.

    idxs = [(0,0,0,0),(0,0,1,1),(0,0,2,2),(0,0,0,1),(0,0,1,2),(0,0,0,2),
            (1,1,0,0),(1,1,1,1),(1,1,2,2),(1,1,0,1),(1,1,1,2),(1,1,0,2),
            (2,2,0,0),(2,2,1,1),(2,2,2,2),(2,2,0,1),(2,2,1,2),(2,2,0,2),
            (0,1,0,0),(0,1,1,1),(0,1,2,2),(0,1,1,0),(0,1,1,2),(0,1,0,2),
            (1,2,0,0),(1,2,1,1),(1,2,2,2),(1,2,0,1),(1,2,2,1),(1,2,0,2),
            (0,2,0,0),(0,2,1,1),(0,2,2,2),(0,2,0,1),(0,2,1,2),(0,2,2,0)]
    """
    M = Matrix([ [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], 
                 [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0] ])
    
    M[0,0], M[0,1], M[0,2] = C[0,0,0,0], C[0,0,1,1], C[0,0,2,2]
    M[0,3], M[0,4], M[0,5] = C[0,0,0,1], C[0,0,1,2], C[0,0,0,2]

    M[1,0], M[1,1], M[1,2] = C[1,1,0,0], C[1,1,1,1], C[1,1,2,2]
    M[1,3], M[1,4], M[1,5] = C[1,1,0,1], C[1,1,1,2], C[1,1,0,2]
    
    M[2,0], M[2,1], M[2,2] = C[2,2,0,0], C[2,2,1,1], C[2,2,2,2]
    M[2,3], M[2,4], M[2,5] = C[2,2,0,1], C[2,2,1,2], C[2,2,0,2]

    #0,1,0,1 -> 0,1,1,0
    M[3,0], M[3,1], M[3,2] = C[0,1,0,0], C[0,1,1,1], C[0,1,2,2]
    M[3,3], M[3,4], M[3,5] = C[0,1,1,0], C[0,1,1,2], C[0,1,0,2]
    
    # 1,2,1,2 -> 1,2,2,1
    M[4,0], M[4,1], M[4,2] = C[1,2,0,0], C[1,2,1,1], C[1,2,2,2]
    M[4,3], M[4,4], M[4,5] = C[1,2,0,1], C[1,2,2,1], C[1,2,0,2]
    
    # 0,2,0,2 -> 0,2,2,0
    M[5,0], M[5,1], M[5,2] = C[0,2,0,0], C[0,2,1,1], C[0,2,2,2]
    M[5,3], M[5,4], M[5,5] = C[0,2,0,1], C[0,2,1,2], C[0,2,2,0]

    return M

def stress_deriv_voigt3D_new(C):
    """
    Given stress tensor derivatives as a fourth order tensor C, 
    this function returns the 3D representation using Voigt notation.

    Reference: Eq.(7.39) J. Bonet    
    """
    M = Matrix([ [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], 
                 [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0] ])
    
    # row 1
    M[0,0] = 2.0 * C[0,0,0,0]
    M[0,1] = 2.0 * C[0,0,1,1]
    M[0,2] = 2.0 * C[0,0,2,2]
    M[0,3] = C[0,0,0,1] + C[0,0,1,0]
    M[0,4] = C[0,0,1,2] + C[0,0,2,1]
    M[0,5] = C[0,0,0,2] + C[0,0,2,0]

    # row 2
    M[1,0] = 2.0 * C[1,1,0,0]
    M[1,1] = 2.0 * C[1,1,1,1]
    M[1,2] = 2.0 * C[1,1,2,2]
    M[1,3] = C[1,1,0,1] + C[1,1,1,0]
    M[1,4] = C[1,1,1,2] + C[1,1,2,1]
    M[1,5] = C[1,1,0,2] + C[1,1,2,0]

    # row 3
    M[2,0] = 2.0 * C[2,2,0,0]
    M[2,1] = 2.0 * C[2,2,1,1]
    M[2,2] = 2.0 * C[2,2,2,2]
    M[2,3] = C[2,2,0,1] + C[2,2,1,0]
    M[2,4] = C[2,2,1,2] + C[2,2,2,1]
    M[2,5] = C[2,2,0,2] + C[2,2,2,0]

    # row 4
    M[3,0] = 2.0 * C[0,1,0,0]
    M[3,1] = 2.0 * C[0,1,1,1]
    M[3,2] = 2.0 * C[0,1,2,2]
    M[3,3] = C[0,1,0,1] + C[0,1,1,0]
    M[3,4] = C[0,1,1,2] + C[0,1,2,1]
    M[3,5] = C[0,1,0,2] + C[0,1,2,0]

    # row 5
    M[4,0] = 2.0 * C[1,2,0,0]
    M[4,1] = 2.0 * C[1,2,1,1]
    M[4,2] = 2.0 * C[1,2,2,2]
    M[4,3] = C[1,2,0,1] + C[1,2,1,0]
    M[4,4] = C[1,2,1,2] + C[1,2,2,1]
    M[4,5] = C[1,2,0,2] + C[1,2,2,0]

    # row 6
    # 3,1 ou 1,3 ?
    M[5,0] = 2.0 * C[2,0,0,0]
    M[5,1] = 2.0 * C[2,0,1,1]
    M[5,2] = 2.0 * C[2,0,2,2]
    M[5,3] = C[2,0,0,1] + C[2,0,1,0]
    M[5,4] = C[2,0,1,2] + C[2,0,2,1]
    M[5,5] = C[2,0,0,2] + C[2,0,2,0]

    return 0.5*M

def print_code_tensor2(S, simplify=True):

    if (simplify):
        detf = Symbol('detf')

    for i in range(3):
        for j in range(3):
            strd = 'S(%d,%d)' % (i,j)
            Sij = S[i,j]

            if (simplify):
                temp = Sij.subs(C.det(),detf*detf)        

            expr = ccode(temp, assign_to=strd)
            print("%s" % expr)


def print_code_tensor4(D, simplify=True):
    nrows, ncols = D.shape

    if (simplify):
        detf = Symbol('detf')

    for i in range(nrows):
        for j in range(ncols):
            strd = 'D(%d,%d)' % (i,j)
            item = D[i,j]

            if (simplify):
                temp = item.subs(C.det(),detf*detf)        

            expr = ccode(temp, assign_to=strd)
            print("%s" % expr)

def print_code_tensor4a(D, simplify=True):
    a,b,c,d = D.shape

    if (simplify):
        detf = Symbol('detf')

    for i in range(a):
        for j in range(b):
            for k in range(c):
                for l in range(d):
                    strd = 'D(%d,%d,%d,%d)' % (i,j,k,l)
                    item = D[i,j,k,l]

                if (simplify):
                    temp = item.subs(C.det(),detf*detf)        

                expr = ccode(temp, assign_to=strd)
                print("%s" % expr)
