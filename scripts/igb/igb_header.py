#!/usr/bin/env python

from numpy import *

"""
IGB file header class

The header is 1024 bytes long and written in ascii. Furthermore the header is divided into lines
of maximum length 80 characters which are separated with new lines characters ('\n'). Each line
contains keyword and value pairs separated from each other by a space. A colon separates the
keyword from its value.

Typical calls for the constructor with arguments:
    
    N = array([400,200,5,1000])
    tipo = 'float'
    dim = array([4,2,0.05,10000])
    uni = ('cm', 'cm', 'cm', 'ms', 'mV')
    hd = igb_header(N,tipo,dim,uni,1,0)
    
    or
    
    hd = igb_header([400,200,5,1000],'float',[4,2,0.05,10000],('cm', 'cm', 'cm', 'ms', 'mV'),1,0)

Bernardo M. Rocha    
"""
class igb_header:
    def __init__(self, N=None, type=None, dim4=None, units=None, facteur=None, zero=None, systeme=None):
        """igb_header constructor --- DEFAULT or with all the ARGUMENTS below"""
        self.x        =  0
        self.y        =  0
        self.z        =  0
        self.t        =  0
        self.dim_x    =  0
        self.dim_y    =  0
        self.dim_z    =  0
        self.dim_t    =  0
        self.unites_x = ''
        self.unites_y = ''
        self.unites_z = ''
        self.unites_t = ''
        self.unites   = ''
        self.type     = 'float'
        self.systeme  = 'little_endian'
        self.facteur  = 1.0
        self.fac_t    = 1.0
        self.zero     = 0
                
        if N != None:
            if len(N) < 4:
                print " ERROR: Input vector N too short (length must be 4)\n"
            self.x  = N[0]
            self.y  = N[1]
            self.z  = N[2]
            self.t  = N[3]
        
        if dim4 != None:
            if len(dim4) < 4:
                print " ERROR: Input vector dimension too short (length must be 4)\n"
            self.dim_x = dim4[0]
            self.dim_y = dim4[1] 
            self.dim_z = dim4[2] 
            self.dim_t = dim4[3]

        if units != None:
            if len(units) < 5:
                print " ERROR: Input vector units too short (length must be 5)\n"
            self.unites_x = units[0]
            self.unites_y = units[1]
            self.unites_z = units[2]
            self.unites_t = units[3]
            self.unites   = units[4]

        if facteur != None: self.facteur = facteur
        if zero    != None: self.zero    = zero
        if systeme != None: self.systeme = systeme

    def __str__(self):
        """display the values of the members of the class"""
        infomsg  = "\t header.x        : "      + str(self.x) + "\n"
        infomsg += "\t header.y        : "      + str(self.y) + "\n"
        infomsg += "\t header.z        : "      + str(self.z) + "\n"
        infomsg += "\t header.t        : "      + str(self.t) + "\n"
        infomsg += "\t header.dim_x    : "    + str(self.dim_x) + "\n"
        infomsg += "\t header.dim_y    : "    + str(self.dim_y) + "\n"
        infomsg += "\t header.dim_z    : "    + str(self.dim_z) + "\n"
        infomsg += "\t header.dim_t    : "    + str(self.dim_t) + "\n"        
        infomsg += "\t header.unites_x : " + str(self.unites_x) + "\n"
        infomsg += "\t header.unites_y : " + str(self.unites_y) + "\n"
        infomsg += "\t header.unites_z : " + str(self.unites_z) + "\n"
        infomsg += "\t header.unites_t : " + str(self.unites_t) + "\n"
        infomsg += "\t header.unites   : "   + str(self.unites) + "\n"
        infomsg += "\t header.type     : "     + str(self.type) + "\n"
        infomsg += "\t header.systeme  : "  + str(self.systeme) + "\n"
        infomsg += "\t header.facteur  : "  + str(self.facteur) + "\n"
        infomsg += "\t header.fac_t    : "    + str(self.fac_t) + "\n"
        infomsg += "\t header.zero     : "     + str(self.zero)        
        return infomsg

    def initialize(self):
        """ initialize header structure --- use NaN for numbers and '' for strings as initializer"""        
        self.x        = NaN
        self.y        = NaN
        self.z        = NaN
        self.t        = NaN
        self.dim_x    = NaN
        self.dim_y    = NaN
        self.dim_z    = NaN
        self.dim_t    = NaN
        self.unites_x = ''
        self.unites_y = ''
        self.unites_z = ''
        self.unites_t = ''
        self.unites   = ''
        self.type     = ''
        self.systeme  = ''
        self.facteur  = NaN
        self.fac_t    = NaN
        self.zero     = NaN
