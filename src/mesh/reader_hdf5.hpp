#ifndef READER_HDF5_HPP
#define READER_HDF5_HPP

#include <string>
#include "mesh.hpp"

/** Reads a mesh from an HDF5 file produced by WriterHDF5.

    Expects the two datasets that WriterHDF5 always writes:
       /geometry/coordinates   (np x 3)  double
       /topology/connectivity  (ne x nn) int

    The element type is inferred from the number of nodes per element
    (nn) together with the spatial dimension, mirroring the mapping used
    in WriterHDF5::write_xdmf.
*/
class ReaderHDF5
{
public:
    //! Read "<file>.h5" (or an explicit .h5 path) and fill mesh.
    //! Set ndim to 2 or 3 so nn=4 can be disambiguated (quad vs tet).
    static void read(const std::string & h5file, Mesh & mesh, int ndim = 3);

    //! Read and create the mesh object
    Mesh * read_mesh(const std::string & h5file, int ndim = 3);
};

#endif // READER_HDF5_HPP
