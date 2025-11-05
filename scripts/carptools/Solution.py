import sys, os
import numpy as np
from scipy import sparse
from scipy.interpolate import *
from fipy import *
from pyUtils.igb.igb_read import *
from pyUtils.petsc.binary_read import petsc_binary_read
from pyUtils.sctools import read_array_pts, read_array_elem

class Solution():
    """
    This class holds the solution data of CARP, such as transmembrane potential (Vm)
    defined on the nodes, FEM matrices (mass and stiffness) and operates on these data.
    """    

    def __init__(self, basename, outdir='output'):

        self.points   = '%s.pts'  % (basename)
        self.elements = '%s.elem' % (basename)

        soldir = os.path.dirname(basename)

        self.vm    = '%s/%s/vm.igb.gz' % (soldir, outdir)
        self.phie  = '%s/%s/phie.igb.gz'% (soldir, outdir)
        self.mass  = '%s/%s/MatLabDump_Mi.bin' % (soldir, outdir)
        self.masse = '%s/%s/MatLabDump_Me.bin' % (soldir, outdir)
        self.massc = '%s/%s/MyMassMatrix_Mi.bin' % (soldir, outdir)
        self.stiff = '%s/%s/MatLabDump_Ki_CN.bin' % (soldir, outdir)

        self.nodes = None

    def getNodes(self):
        """
        Reads in the points file and return an array of the form a(numPoints,3)
        """
        if self.nodes is None:
            if os.path.isfile(self.points):
                self.nodes = read_array_pts(self.points)
                return self.nodes
        else:
            return self.nodes    

    def getNumberOfNodes(self):
        """
        Returns the number of nodes of the solution's mesh
        """
        temp = self.getNodes()
        return len(temp)
        
    def getElements(self):
        """
        Reads in the element file and return an array/list of elements
        """
        if os.path.isfile(self.elements):
            return read_array_elem(self.elements)

    def getVm(self):
        """
        Reads in IGB file and returns Vm array
        """
        if os.path.isfile(self.vm):          
            [vm, header] = read_igb_slice(self.vm, is_gzipped=True)
            return vm.squeeze()   

    def getPhie(self):
        """
        Reads in IGB file and returns Vm array
        """
        if os.path.isfile(self.phie):          
            [phie, header] = read_igb_slice(self.phie, is_gzipped=True)
            return phie.squeeze()

    def getLumpedMassMatrix (self):
        """
        Reads in PETSc binary vector for the lumped mass matrix and returns it as
        a scipy.sparse COO matrix
        """
        if not os.path.isfile(self.mass):
            print ' Error: the file %s does not exist.' % self.mass
            sys.exit(-1)
        
        binaryfile = self.mass
        data = petsc_binary_read (binaryfile,0)
        size = np.size(data)
        rows = np.arange(0,size,1)
        cols = np.arange(0,size,1)
        mass = sparse.coo_matrix( (data,(rows,cols)) ,(size, size) )
        return mass

    def getLumpedMassMatrixE (self):
        """
        Reads in PETSc binary vector for the lumped mass matrix and returns it as
        a scipy.sparse COO matrix
        """
        if not os.path.isfile(self.masse):
            print ' Error: the file %s does not exist.' % self.masse
            sys.exit(-1)
        
        binaryfile = self.masse
        data = petsc_binary_read (binaryfile,0)
        size = np.size(data)
        rows = np.arange(0,size,1)
        cols = np.arange(0,size,1)
        mass_e = sparse.coo_matrix( (data,(rows,cols)) ,(size, size) )
        return mass_e

    def getMassMatrix (self):
        """
        Reads in PETSc binary matrix for the mass matrix and returns it as a scipy.sparse COO matrix
        """
        if not os.path.isfile(self.mass):
            print ' Error: the file %s does not exist.' % self.mass
            sys.exit(-1)
            
        binaryfile = self.mass
        return petsc_binary_read(binaryfile,0)

    def getStiffnessMatrix(self): 
        """
        Reads in PETSC binary matrix for the stiffness matrix
        """
        if not os.path.isfile(self.stiff):
            print ' Error: the file %s does not exist.' % self.stiff
            sys.exit(-1)

        binaryfile = self.stiff
        return petsc_binary_read(binaryfile,0)

    def getL2Norm (self):
        """
        Compute the L2 Norm of two (Vm) IGB files
        """
        # reads Vm
        vm1 = getVm (igbfile1)
        vm2 = getVm (igbfile2)
        
        # reads CARP lumped mass matrix and store in scipy.sparse COO matrix
        mass1 = getMassMatrix (massmatrix1)
        mass2 = getMassMatrix (massmatrix2)
        
        # compute norm
        nrm1 = getL2NormError (vm1[time,:],mass1)
        nrm2 = getL2NormError (vm2[time,:],mass2)
        
        return abs(nrm1-nrm2)

    def getL2NormError (self, ref, time):
        """
        Compute the L2 norm of two (Vm) IGB files given the reference solution.
        The coarse solution is interpolated to the reference mesh and the the L2 norm is computed 
        """
        # reads in Vm
        vmsol = self.getVm()
        vmref = ref.getVm()

        # reads in mass matrix
        mass1 = ref.getLumpedMassMatrix()
        mass2 = self.getLumpedMassMatrix()

        # read mesh
        meshsol = self.getNodes()
        meshref = ref.getNodes()
        
        vm1 = vmref[time,:]
        vm2 = vmsol[time,:]
        
        # interpolate to reference mesh
        vm2i = self.intepolateVmToRef(meshref, meshsol, vm2)
        
        # compute error as a column
        e = vm1[:] - vm2i[:]
        
        # compute the relative L2 norm
        norm0 = self.calcL2NormError(e,mass1) / self.calcL2NormError(vm1,mass1)
        return norm0

    def getMaxNormError(self, ref, time):
        """
        Compute the Max norm of two (Vm) IGB files interpolating the solution to the 
        reference mesh
        """
        # read in Vm
        vmref = ref.getVm()
        vmsol = self.getVm()

        # read in mass matrix
        mass1 = ref.getLumpedMassMatrix()
        mass2 = self.getLumpedMassMatrix()

        # read mesh
        meshsol = self.getNodes()
        meshref = ref.getNodes()

        vm1 = vmref[time,:]
        vm2 = vmsol[time,:]

        # interpolate to reference mesh
        vm2i = self.intepolateVmToRef(meshref, meshsol, vm2)

        # compute error
        e = np.abs(vm1[:] - vm2i[:])
        return np.max(e)

    def calcL2NormError(self, e, mass):
        """
        Compute the L2 norm of the error vector e, given the Mass matrix mass
        """
        aux = np.dot(e, mass*e)
        return np.sqrt(aux)
        
    def getArrivalTime (self):
        """
        Compute the arrival time for one node using the maximum dVm/dt
        """
        vm = self.getVm()
        return np.max( np.diff(vm, n=1) )

    def intepolateVmToRef(self, meshref, meshsol, vmsol):
        """
        Interpolate Vm from a coarse mesh to a fine mesh using FiPy routines
        Take into account that it is a structured mesh generated by mesher.

        Depends on package ***FiPY***
        """
        # define fine mesh\textbf{
        dx = meshref[1,0] - meshref[0,0]
        dy = dx
        nx = int( ((meshref[-1,0]*2)/dx ) + 1 )
        ny = int( ((meshref[-1,1]*2)/dy ) + 1 )    
        mesh1 = Grid2D (nx=nx, ny=ny, dx=dx, dy=dy)
        
        # define coarse mesh
        dx = meshsol[1,0] - meshsol[0,0]
        dy = dx
        nx = int( (( abs(meshsol[0,0])*2)/dx) + 1 )
        ny = int( (( abs(meshsol[0,1])*2)/dy) + 1 )    
        mesh2 = Grid2D (nx=nx, ny=ny, dx=dx, dy=dy)
        
        # define variable vm over coarse grid
        vm = CellVariable(mesh=mesh2, value=vmsol)
        
        # interpolate vm from coarse to fine grid
        vm_sol_int = np.array(vm(mesh1.getCellCenters(), order=1))
        
        return vm_sol_int
#
# end of class Solution
#
