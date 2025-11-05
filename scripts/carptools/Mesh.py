#!/usr/bin/env python

import numpy as np
from pyUtils.sctools import read_array_pts, read_array_elem

elemInfo = []

elemInfo.append( {'label': 'Tt', 'numNodes': 4, 'numEdges': 6,  'numFaces': 4, 'dim': 3} )
elemInfo.append( {'label': 'Hx', 'numNodes': 8, 'numEdges': 12, 'numFaces': 6, 'dim': 3} )
elemInfo.append( {'label': 'Oc', 'numNodes': 6, 'numEdges': 12, 'numFaces': 8, 'dim': 3} )
elemInfo.append( {'label': 'Py', 'numNodes': 5, 'numEdges': 8,  'numFaces': 5, 'dim': 3} )
elemInfo.append( {'label': 'Pr', 'numNodes': 6, 'numEdges': 9,  'numFaces': 5, 'dim': 3} )
elemInfo.append( {'label': 'Qd', 'numNodes': 4, 'numEdges': 4,  'numFaces': 1, 'dim': 2} )
elemInfo.append( {'label': 'Tr', 'numNodes': 3, 'numEdges': 3,  'numFaces': 1, 'dim': 2} )
elemInfo.append( {'label': 'Ln', 'numNodes': 2, 'numEdges': 2,  'numFaces': 0, 'dim': 1} )

class CarpMesh:
    """
    Class for carp meshes
    """

    def __init__(self, basename):

        self.baseFile     = basename
        self.pointsFile   = '%s.pts'  % basename
        self.elementsFile = '%s.elem' % basename
    
        pts  = read_array_pts (self.pointsFile)
        elem = read_array_elem (self.elementsFile)

        self.points   = pts
        self.elements = elem

        self.numPoints   = np.shape(self.points)[0]
        self.numElements = np.shape(self.elements)[0]    

    def getNumPoints (self):
        """
        Returns the number of points/nodes in the mesh
        """
        return self.numPoints
    
    def getNumElements (self):
        """
        Returns the number of elements in the mesh
        """
        return self.numElements

    def getNode (self, i):
        """
        Returns the points of the node i
        """
        return self.points[n]

    def getElementNodes (self, i):
        """
        Returns the element global nodes
        """
        elemType = self.elements[i][0]
        print elemType

   
if __name__ == "__main__":
    
    m = CarpMesh ('/home/rocha/src/pyUtils/testdata/tri_100um')
    
    print m.getNumPoints()
    print m.getNumElements()
    
    m.getElementNodes(500)

