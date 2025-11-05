#!/usr/bin/env python

from pyUtils.carptools.ParFile import ParameterFile

class Cable(ParameterFile):
    def __init__(self, im, gil, vs):
        ParameterFile.__init__(self, ionicModel=im, carp_ver=vs)

        self.set_parameter('solnmethod',  4)
        self.set_parameter('gridout_i',   2)
        self.set_parameter('dt',          10)
        self.set_parameter('timedt',      0.5)
        self.set_parameter('spacedt',     1)
        self.set_parameter('CN_parab',    1)
        self.set_parameter('cg_tol_parab',1.0e-06)
        self.set_parameter('tend',        200)

        # version specific settings
        if vs is 'carpe':
          self.set_parameter('readmesh',  3)
          self.set_parameter('gil',       gil)

        if vs is 'carpm':
          self.set_parameter('gregion[0].g_il', gil)

        self.add_LAT ('activation')

# end of Cable

class Cable2D(Cable):
    def __init__(self, im, gil, vs):
        ParameterFile.__init__(self, ionicModel=im, carp_ver=vs)
        
        self.add_stimulus (0, 0, 5e-2, 1, 100, 1000, 100, -10099, -500, 0)
        
# end of Cable2D

class Cable3D(Cable):
    def __init__(self, im, gil, vs):
        ParameterFile.__init__(self, ionicModel=im, carp_ver=vs)

        self.add_stimulus (0, 0, 5e-2, 1, 100, 1000, 1000, -10099, -500, -500)
        
# end of Cable3D