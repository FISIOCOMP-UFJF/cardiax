import sys, inspect
import numpy as np
from sympy import Symbol
from mechutils import *
from hyperelastic_base import *

class StVenantKirchoff(HyperelasticMaterial):
    def __init__(self):
        self.mu = Symbol('mu')
        self.lb = Symbol('lambda')

    def strain_energy(self):
        mu, lb = self.mu, self.lb
        E = 0.5*(C - I) 
        return 0.5*lb*(trace(E)**2) + mu*contract(E,E)       

class CompressibleNeoHookean(HyperelasticMaterial):

    def __init__(self):
        self.mu = Symbol('mu')
        self.lb = Symbol('lambda')

    def strain_energy(self):
        mu, lb = self.mu, self.lb
        detF = sqrt(I3(C))
        return (mu/2.0)*(I1(C) - 3.0) - mu * log(detF) + \
            (lb/2.0)*(log(detF)*log(detF))

class CompressibleMooneyRivlin(HyperelasticMaterial):
    
    def __init__(self):
        self.c1 = Symbol('c1')
        self.c2 = Symbol('c2')
        self.K  = Symbol('K')        

    def strain_energy(self):
        c1, c2, K = self.c1, self.c2, self.K
        detF    = sqrt(I3(C))
        I1C_hat = I1_deviatoric(C)
        I2C_hat = I2_deviatoric(C)
        return c1*(I1C_hat-3.0) + c2*(I2C_hat-3.0) + (K/2.0)*(log(detF)**2)

class SCompressibleMooneyRivlin(HyperelasticMaterial):
    
    def __init__(self):
        self.c1 = Symbol('c1')
        self.c2 = Symbol('c2')
        self.Ccompr = Symbol('C_Compr')

    def strain_energy(self):
        c1, c2, Ccompr = self.c1, self.c2, self.Ccompr
        I3C = sqrt(I3(C)) # detF
        I1C = I1(C)
        I2C = I2(C)
        return c1*(I1C-3.0) + c2*(I2C-3.0) -(c1+2.0*c2)*log(I3C) + Ccompr*(I3C*log(I3C)-I3C+1.0)

class FungType(HyperelasticMaterial):

    def __init__(self):
        self.K = Symbol('K')
        self.Ccompr = Symbol('Ccompr')
        self.bff = Symbol('bff')
        self.bxx = Symbol('bxx')
        self.bfx = Symbol('bfx')

    def strain_energy(self):
        f,s,n = 0, 1, 2
        K,bff,bxx,bfx = self.K, self.bff, self.bxx, self.bfx
        Ccompr = self.Ccompr
        J = sqrt(I3(C))
        E = 0.5 * (C - I)
        Q  = bff * (E[f,f]**2)
        Q += bxx * (E[n,n]**2 + E[s,s]**2 + E[s,n]**2 + E[n,s]**2)
        Q += bfx * (E[f,n]**2 + E[n,f]**2 + E[f,s]**2 + E[s,f]**2)
        return 0.5 * K * (exp(Q) - 1.0) + Ccompr*(J*log(J)-J+1)

class Usyk(HyperelasticMaterial):
    
    def __init__(self):
        self.Kfu = Symbol('Kfu')
        self.bff = Symbol('bff')
        self.bxx = Symbol('bxx')
        self.bfx = Symbol('bfx')
        self.c_param = Symbol('c_param')
        self.c_compr = Symbol('c_compr')

    def strain_energy(self):
        E = 0.5 * (C - I)
        fw = bff*E[0,0]*E[0,0]
        fw  += bxx*(E[1,1]*E[1,1] + E[2,2]*E[2,2] + E[1,2]*E[1,2] + E[2,1]*E[2,1])
        fw  += bfx*(E[0,2]*E[0,2] + E[2,0]*E[2,0] + E[0,1]*E[0,1] + E[1,0]*E[1,0])
        psi  = 0.5 * c_param * (exp(fw - 1.0)) 
        psi += c_compr * (J * log(J) - J + 1.0)
        return psi

class PoleZero(HyperelasticMaterial):
    
    def __init__(self):
        self.k11, self.k22, self.k33 = symbols('k11,k22,k33')
        self.k12, self.k13, self.k23 = symbols('k12,k13,k23')
        self.a11, self.a22, self.a33 = symbols('a11,a22,a33')
        self.a12, self.a13, self.a23 = symbols('a12,a13,a23')

    def strain_energy(self):
        t1 = k11 * ((e11*e11)/(Abs(a11-e11)))
        t2 = k22 * ((e22*e22)/(Abs(a22-e22)))
        t3 = k33 * ((e33*e33)/(Abs(a33-e33)))
        t4 = k12 * ((e12*e12)/(Abs(a12-e12)))
        t5 = k13 * ((e13*e13)/(Abs(a13-e13)))
        t6 = k23 * ((e23*e23)/(Abs(a23-e23)))
        W = t1 + t2 + t3 + t4 + t5 + t6
        return W

def create_hyperelastic_material(name):
    """
    Factory function to create differente Hyperelastic materials
    """
    if   (name == 'CompressibleNeoHookean'): return CompressibleNeoHookean()
    elif (name == 'CompressibleMooneyRivlin'): return CompressibleMooneyRivlin()
    elif (name == 'SCompressibleMooneyRivlin'): return SCompressibleMooneyRivlin()
    elif (name == 'StVenantKirchoff'): return StVenantKirchoff()
    elif (name == 'FungType'): return FungType()
    else:
        print("Error: unknown material type.")
        sys.exit(1)

