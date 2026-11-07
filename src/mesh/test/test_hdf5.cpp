#include <cstdlib>
#include "../mesh.hpp"
#include "../writer_hdf5.hpp"

using namespace std;

int main(int argc, char** argv)
{
    cout << "Testing HDF5 + XDMF writer" << endl;

    int nsteps = 20;
    std::string nome = argv[1];
    Mesh msh;
    msh.read(nome);

    WriterHDF5 out(&msh);

    // register fields BEFORE open() so the file holds only these
    out.add_field("vm",            FieldLoc::Point, FieldKind::Transient);
    out.add_field("displacements", FieldLoc::Point, FieldKind::Transient, 3);
    out.open("output_test_xdmf", nsteps, 1.0/nsteps);

    double *v = new double[msh.get_n_points()];
    double *u = new double[3*msh.get_n_points()];

    for(int i=0; i<nsteps; i++)
    {
        for(uint j=0; j<msh.get_n_points(); j++)
        {
            v[j] = i*10;                 // scalar field
            u[j*3 + 0] = i*200;          // displacement field
            u[j*3 + 1] = i*500;
            u[j*3 + 2] = i*1000;
        }

        out.write_point_field_step(i, v, "vm");   // was write_vm_step(i, v)
        out.write_displ_step(i, u);               // unchanged
    }

    delete [] v;
    delete [] u;                         // <-- was leaked; added

    out.close();

    cout << "Done" << endl;
    return 0;
}