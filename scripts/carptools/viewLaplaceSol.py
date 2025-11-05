#!/usr/bin/env python

import os, sys
import numpy as np
from pylab import *
from pyUtils.carptools.Solution import Solution

def viewLaplaceSolution (solname):

    lasol = Solution(solname)
    lphie = lasol.getPhie()

    # no time steps
    lphie = lphie[0,:]
    
    peshp = np.shape(lphie)[0]

    # view structured square mesh only
    numNodesX = int( sqrt(peshp) )
    numNodesY = int( sqrt(peshp) )

    # reshape for visualization
    phie = np.reshape(lphie, (numNodesX,numNodesY))

    ax = subplot(111)
    im = imshow(phie,cmap=cm.jet)
    im.set_interpolation('bilinear')
    title('Solution of the Laplace problem')
    colorbar()
    show()

if __name__ == "__main__":

    if len(sys.argv) > 1:
        viewLaplaceSolution(sys.argv[1])
    else:
        print " Usage: viewLaplaceSol <SolutionDir> \n"
    
