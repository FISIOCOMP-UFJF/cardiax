import sys

if (len(sys.argv) < 3):
    print(" Usage: python convelastic.py <E> <nu>")
    exit(1)

# Young modulus
E  = float(sys.argv[1])

# Poisson ratio
nu = float(sys.argv[2])

# Lame parameters
mu   = E/(2.0*(1.0+nu))
lamb = (E*nu)/((1 + nu)*(1 - 2*nu))
K = (E/(3.0*(1-2.0*nu)))

print("\nElasticity parameters")
print(" Young modulus E: %f" % E)
print(" Poisson ratio Nu: %f" % nu)
print("\nLame parameters")
print("Mu            Lambda        Bulk modulus")
print("%.6e  %.6e  %.6e" % (mu,lamb,K))



