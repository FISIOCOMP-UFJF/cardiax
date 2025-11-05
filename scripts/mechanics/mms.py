import sys, inspect
from sympy import Matrix
from hyperelastic_materials import *

def usage():
    ustr = """
 Usage: mms <material>
 
  Options:
    material: name of the material (NeoHookean, MooneyRivlin, StVenantKirchoff)
"""
    print(ustr)
# end of usage
    
def PiolaKirchoff2(E):
    lmbda, mu = symbols('lambda mu')
    return lmbda*E.trace()*I + 2*mu*E # St. Venant-Kirchoff model

def WriteCode(f):
    print("*** BEGIN ***")
    strcode = '''
  Vector f(3);
  Real x = X(0), y = X(1), z = X(2);
  Real k = 2.0;
  Real mu = mat.get_param(0);
  Real lambda = mat.get_param(1);
  f(0) = %(f0)s;
  f(1) = %(f1)s;
  f(2) = %(f2)s;
  return(f);
''' % dict(f0=ccode(f[0]),f1=ccode(f[1]),f2=ccode(f[2]))
    print(strcode)
    print("*** END ***")

if __name__ == "__main__":

    x,y,z = symbols('x y z')
    X = Matrix([x,y,z]) 
    
    # Solution
    k = Symbol('k')
    #u = Matrix([ k*sin(x), 0,  0 ])
    #u = Matrix([ k*y*y, 0,  0 ])

    u = Matrix([x**4 + 2**y*z/5, y**4 + 2*x*z/5, z**4/10 - 2*x*y*z/5])

    # Kinematics
    du = u.jacobian(X)
    H = du
    F = I + H
    C = F.T * F
    E = 0.5 * (C - I)
    S = PiolaKirchoff2(E)
    Pi = F * S

    # Compute strong-form
    forcingvec = [-divergence(row,X) for row in rows(Pi)]

    # Print data
    print("\nGrad u")
    print(H)
    print("\nDeformation gradient F")
    print(F)
    print("\nRight Cauchy-Green strain tensor")
    print(C)
    print("\nLagrangian strain tensor E")
    print(E)
    print("\nPiola-Kirchoff 2")
    print(S)
    print("\nPiola-Kirchoff 1")
    print(Pi)

    print("\ndiv P")
    print(forcingvec)      
    
    WriteCode(forcingvec)

    print("\nDone")
