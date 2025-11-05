# Some utils
import numpy as np
from math import sin, cos, tan

def rotateAroundAxis(v,u,theta):
    """
    Rotate v by an angle of theta about an axis in the direction of u.
    """
    #theta = theta*(pi/180)

    R = np.zeros((3,3))
    c = cos(theta)
    s = sin(theta)
    ux,uy,uz = u[0], u[1], u[2]
    ux2 = ux*ux
    uy2 = uy*uy
    uz2 = uz*uz

    R[0,0] = c + ux2*(1.-c)
    R[0,1] = ux*uy*(1.-c) - uz*s
    R[0,2] = ux*uz*(1.-c) + uy*s

    R[1,0] = uy*ux*(1.-c) + uz*s
    R[1,1] = c + uy2*(1.-c)
    R[1,2] = uy*uz*(1.-c) - ux*s

    R[2,0] = uz*ux*(1.-c) - uy*s
    R[2,1] = uz*uy*(1.-c) + ux*s
    R[2,2] = c + uz2*(1.-c)

    R = np.matrix(R)
    v = np.matrix([[v[0]],[v[1]],[v[2]]])

    return np.array(R * v).reshape(3,)   
