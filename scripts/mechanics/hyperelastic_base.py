import sys
from sympy import *
from sympy.printing import print_ccode
from sympy.abc import x, y, z
import numpy as np

# Right Cauchy-Green strain tensor
c11, c12, c13 = symbols('c11,c12,c13')
c21, c22, c23 = symbols('c21,c22,c23')
c31, c32, c33 = symbols('c31,c32,c33')

C = Matrix([ [c11,c12,c13], [c21,c22,c23], [c31,c32,c33] ])
I = Matrix([ [1,0,0], [0,1,0], [0,0,1] ])
#I = eye(3)

class HyperelasticMaterial():
    """
    Base class for hyperelastic materials
    """

    def __init__(self):
        self.num_parameters = 0
        self.parameters = []

    def strain_energy(self):
        pass

    def calc_stress(self, type="Compressible"):
        W = self.strain_energy
        S = Matrix([ [0,0,0], [0,0,0], [0,0,0] ])    

        if (type == "Compressible"):
            for i in range(3):
                for j in range(3):
                    S[i,j] = 2.0 * diff(W(), C[i,j])
            return S

        elif (type == "NearlyIncompressible"):
            pass

    def calc_stress_deriv(self, type="Compressible"):
        W = self.strain_energy
        S_deriv = np.zeros((3,3,3,3),dtype=object)
        
        if (type == "Compressible"):
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            S_deriv[i,j,k,l] = 4.0 * diff(W(), C[i,j], C[k,l])
            return S_deriv
    

    
