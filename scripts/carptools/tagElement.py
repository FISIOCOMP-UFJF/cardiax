#!/usr/bin/env python

import sys, gzip, pdb
from numpy import loadtxt, argmax, argmin, array, where, zeros
from sctools import read_array_elem, read_array_pts

class Point:
    
    def __init__(self, x, y, z):
        self.x = x; self.y = y; self.z = z;
    
    def setX(self,X):
        self.x = X

    def setY(self,Y):
        self.y = Y

    def setY(self,Z):
        self.z = Z
    
    def __str__(self):
        return 'x=%8.3f y=%8.3f z=%8.3f' % (self.x, self.y, self.z)

# end of class Point

class RectangularRegion():

    def __init__(self, b, t, w, h, tag):
        self.bottom = b
        self.top    = t
        self.width  = w
        self.height = h
        self.tag    = tag
    
    def __str__(self):
        retstr  = 'left bottom: %f %f' % (self.bottom.x, self.bottom.y)
        retstr += 'right top  : %f %f' % (self.top.x, self.top.y)
        retstr += 'width      : %f' % self.width
        retstr += 'height     : %f' % self.height
        return retstr

# end of class RectangularRegion

def tag_elem_file (elem_file, node_file, region_list):
    """
    Function to tag rectangular areas in structured meshes in CARP
    Bernardo M. Rocha, 2009
    """   
    new_efile = open('element.elem', 'w')
    
    elem = read_array_elem (elem_file)
    pts  = read_array_pts  (node_file)

    for r in region_list:
        
        for e in elem:
            elem_type = e[0]   
            if elem_type == 'Qd':
                node_range = xrange(1,4+1)
                node_tag   = e[5]
                for n in node_range:
                    g = int( e[n] )
                    node = Point ( pts[g][0], pts[g][1], pts[g][2] )        
                    # check if node is inside region
                    if ( (node.x >= r.bottom.x) and
                         (node.x <= r.bottom.x + r.width) and
                         (node.y >= r.top.y - r.height) and
                         (node.y <= r.top.y)):
                        node_tag = r.tag
                    # end if
                # end for n
                e[5] = node_tag
            # end if Qd
        # end for e
    # end for r
    
    # write to file
    
    new_efile.write('%d\n' % len(elem))
    for e in elem:
        new_efile.write('Qd %d %d %d %d %d \n' %
            (int(e[1]), int(e[2]),int(e[3]),int(e[4]),int(e[5])))

    new_efile.close()

    new_efile.close()
    
    

if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print "\n Usage: carp_tagelem.py [elem_list] [node_list]\n"    
        exit(-1)
    
    # regions width and height in um
    w = 2000.0; h = 2000.0
    
    # region 1 in mm
    pt1  = Point(-3000.0, -1000.0, 0.0) # left bottom
    pt2  = Point(-1000.0, -1000.0, 0.0) # right top
    reg1 = RectangularRegion(pt1, pt2, w, h, 10)
    
    # region 2 in mm
    pt1  = Point( 1000.0,  3000.0, 0.0) # left bottom
    pt2  = Point( 3000.0, -1000.0, 0.0) # right top
    reg2 = RectangularRegion(pt1, pt2, w, h ,20)
    
    # tag a list of regions
    region_list = []  
    region_list.append(reg1)
    region_list.append(reg2)
    
    tag_elem_file(sys.argv[1], sys.argv[2], region_list)
                
    #pdb.set_trace()


