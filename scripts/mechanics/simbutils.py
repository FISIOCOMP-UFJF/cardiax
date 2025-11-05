import sys
import numpy as np
import sympy

def grad(U, X):
    return U.jacobian(X)

def symgrad(U, X):
    return sym(grad(U,X))

def sym(A):
    return 0.5*(A + A.T)

def divergence(F, X):
    assert(F.cols == X.rows)
    return sum(sympy.diff(F[i],X[i]) for i in range(X.rows))

def symbols(vars):
    sym = sympy.symbols(list(vars))
    try: # make sure that result is a list even if it only contains one entry
        sym[0]
        return sym
    except:
        return [sym]

def symbol3(var):
    return symbols(['%s[%d]' % (var, i) for i in range(3)])

def symbol33(var):
    return [symbol3('%s[%d]' % (var,i)) for i in range(3)]

if __name__ == "__main__":
    pass
