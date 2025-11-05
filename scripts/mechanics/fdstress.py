# -*- coding: utf-8 -*-

from numpy import *

eps = 1.0e-4

I = matrix([[1.,0,0],[0,1.,0],[0,0,1.]])

def trace(T):
    return T[0,0]+T[1,1]+T[2,2]

def contract(A,B):
    s = 0
    for i in range(3):
        for j in range(3):
            s += A[i,j] * B[i,j]
    return s

def kronecker(i,j):
    return (1 if i==j else 0)

def deltaTensor(k,l):
    """
    Computes the symmetric pertubation 'delta tensor' for the strain energy
    """
    D = zeros((3,3))
    for i in range(3):
        for j in range(3):
            D[i,j] = 0.5*(kronecker(i,k)*kronecker(j,l)+
                           kronecker(i,l)*kronecker(j,k))*eps        
    return D

def strain_energy(E):
    mu, lmbda = 0.3846, 0.8333
    return 0.5*lmbda*(trace(E)**2) + mu*contract(E,E)

def svk_stress(mu,lambd,C):
    E = 0.5*(C-I)
    return lambd*trace(E)*I + 2*mu*E

def svk_elastens(mu,lambd):
    C = zeros((3,3,3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3): 
                    C[i,j,k,l] = lambd*kronecker(i,j)*kronecker(k,l) 
                    C[i,j,k,l] += mu*kronecker(i,k)*kronecker(j,l)
                    C[i,j,k,l] += mu*kronecker(i,l)*kronecker(j,k)
    return C

if __name__ == "__main__":
    
    S = zeros((3,3))
    C = zeros((3,3))
    C[0,0] = 2.
    C[0,1] = 1
    C[1,1] = 2
    C[1,0] = 1
    C[2,2] = 2

    E = 0.5*(C-identity(3))
    
    for i in range(3):
        for j in range(3):
            delta = deltaTensor(i,j)
            psip  = strain_energy(E + delta)
            psim  = strain_energy(E - delta)
            S[i,j] = (psip - psim)/(2.0*eps)

    print("Finite Difference - PK2")
    print(S)

    print("Analytical - PK2")
    mu, lmbda = 0.3846, 0.8333
    print(svk_stress(mu,lmbda,C))

    D = zeros((3,3,3,3))
    for i in range(3):
        for j in range(3):
            delta_ij = deltaTensor(i,j)
            for k in range(3):
                for l in range(3):    
                    delta_kl = deltaTensor(k,l)

                    if(i==k and j==l):
                        psip = strain_energy(E + delta_ij)
                        psim = strain_energy(E - delta_ij)
                        psi0 = strain_energy(E) 
                        D[i,j,k,l] = (psip - 2*psi0 + psim)/(eps**2)
                        print(psip,psim,psi0,eps**2)
                    else:
                        psipp = strain_energy(E + delta_ij + delta_kl)
                        psimp = strain_energy(E - delta_ij + delta_kl)
                        psipm = strain_energy(E + delta_ij - delta_kl)
                        psimm = strain_energy(E - delta_ij - delta_kl)
                        D[i,j,k,l] = (psipp - (psimp+psipm) + psimm)/(4.0*eps**2)

    print("Finite Difference - Material elasticity tensor")
    print(D)

    print("Analytical - Material elasticity")
    M = svk_elastens(mu,lmbda)
    print(M)

    #print("Diff ")
    #print(elastens(mu,lmbda) - D)
    
    #for i in range(3):
    #    for j in range(3):
    #        for k in range(3):
    #            for l in range(3):   
    #                print("%f    %f = %f" % (M[i,j,k,l], D[i,j,k,l], M[i,j,k,l]-D[i,j,k,l]))
