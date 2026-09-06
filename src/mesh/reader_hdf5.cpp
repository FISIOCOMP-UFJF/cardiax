#include "reader_hdf5.hpp"
#include <vector>
#include <stdexcept>
#include <highfive/highfive.hpp>

void ReaderHDF5::read(const std::string & h5file, Mesh & mesh, int ndim)
{
    // Accept either "name" or "name.h5"
    std::string path = h5file;
    if (path.size() < 3 || path.substr(path.size() - 3) != ".h5")
        path += ".h5";

    HighFive::File file(path, HighFive::File::ReadOnly);

    // geometry: (np x 3) doubles
    std::vector<std::vector<double>> coords;
    file.getDataSet("/geometry/coordinates").read(coords);

    const std::size_t np = coords.size();

    // topology: (ne x nn) ints
    std::vector<std::vector<int>> connec;
    file.getDataSet("/topology/connectivity").read(connec);

    const std::size_t ne = connec.size();
    const int nn = ne ? static_cast<int>(connec[0].size()) : 0;

    // element type from (nn, ndim), same logic as write_xdmf
    ElementType etype;
    if      (nn == 2)               etype = ELEM_SEGM;
    else if (nn == 3)               etype = ELEM_TRIG;
    else if (nn == 4 && ndim == 2)  etype = ELEM_QUAD;
    else if (nn == 4 && ndim == 3)  etype = ELEM_TETRA;
    else if (nn == 8)               etype = ELEM_HEXA;
    else throw std::runtime_error("ReaderHDF5: unsupported element (nn=" +
                                  std::to_string(nn) + ")");

    // populate the mesh
    mesh.set_ndim(ndim);
    mesh.set_nen(nn);
    mesh.set_npoints(np);
    mesh.set_nel(ne);
    mesh.reserve_points(np);
    mesh.reserve_elements(ne);

    for (std::size_t i = 0; i < np; ++i)
        mesh.add_point(arma::vec3{coords[i][0], coords[i][1], coords[i][2]});

    // marker/aha default to 0 (not stored in the base mesh datasets).
    for (std::size_t e = 0; e < ne; ++e)
        mesh.add_elem(Element(etype, connec[e], 0, 0));
}

Mesh * ReaderHDF5::read_mesh(const std::string & h5file, int ndim)
{
    Mesh * mesh = new Mesh();
    read(h5file, *mesh, ndim);   // reuse the fill-in-place version
    return mesh;
}