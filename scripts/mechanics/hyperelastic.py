import sys, inspect
from simbutils import *
from hyperelastic_materials import *

def usage():
    ustr = """
 Usage: hyperelastic <material> <dim>
 
  Options:
    material:    name of the material
    dim:         number of dimensions

  Materials:
    FungType
    StVenantKirchoff
    CompressibleNeoHookean
    CompressibleMooneyRivlin
    SCompressibleMooneyRivlin (Slightly compressible - Ref: Usyk)
"""
    print(ustr)
# end of usage

if __name__ == "__main__":

    if (len(sys.argv) < 3):
        usage()
        sys.exit(1)

    mat = sys.argv[1]
    dim = sys.argv[2]

    m   = create_hyperelastic_material(mat)
    S   = m.calc_stress()
    dS  = m.calc_stress_deriv()

    if (dim =='2' or dim == '2D'):
        dSv = stress_deriv_voigt2D(dS)
        dSv2 = stress_deriv_voigt2D_new(dS)
    elif (dim == '3' or dim == '3D'):
        dSv = stress_deriv_voigt3D(dS)
        dSv2 = stress_deriv_voigt3D_new(dS)
    
    print("\nSecond Piola-Kirchhoff stress tensor\n")
    print_code_tensor2(S)

    print("\nFourth order elasticity tensor\n")
    print_code_tensor4(dSv2)

    print("\nDone")
