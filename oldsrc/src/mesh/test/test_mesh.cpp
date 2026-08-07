#include <iostream>

#include "../mesh.hpp"
#include "../box_mesh.hpp"
#include "../purkinje_mesh.hpp"
#include "../writer_hdf5.hpp"

int main()
{
  // --------------------------------------------------------------------------
  // Test 1
  // --------------------------------------------------------------------------

  cout << "TEST 1" << endl;
  /*
  Mesh mesh;
  mesh.read("/data/meshes/heart/benchmark/slab100um");

  std::vector<int> ptnums;
  mesh.get_element_pt_nums(0, ptnums);

  cout << "Elemento 1" << endl;
  cout << ptnums[0] << "\t" << ptnums[1] << "\t"
       << ptnums[2] << "\t" << ptnums[3] << endl;

  cout << "Nos do elemento 1" << endl;
  cout << mesh.get_point(ptnums[0]).t();
  cout << mesh.get_point(ptnums[1]).t();
  cout << mesh.get_point(ptnums[2]).t();
  cout << mesh.get_point(ptnums[3]).t();
  */

  // --------------------------------------------------------------------------
  // Test 2: create a (unit) box mesh
  // --------------------------------------------------------------------------

  cout << endl;
  cout << "TEST 2" << endl;
  cout << "Creating bounding box mesh" << endl;
  BoxMesh bb(0,10,0,10,0,10,2);
  cout << bb << endl;

  cout << "Writing bounding box mesh" << endl;
  int nsteps = 2;

  WriterHDF5 writer(&bb);
  writer.open("testeMesh", nsteps, 1.0/nsteps);

  double *v = new double[bb.get_n_points()];
  double *u = new double[3*bb.get_n_points()];

  for(int i=0; i<nsteps; i++)
  {
    for(uint j=0; j<bb.get_n_points(); j++)
    {
      // scalar field
      v[j] = i*10;

      // displacement field
      u[j*3 + 0] = i*200;
      u[j*3 + 1] = i*500;
      u[j*3 + 2] = i*1000;
    }

    writer.write_vm_step(i, v);
    writer.write_displ_step(i, u);
  }

  delete [] v;
  writer.close();

  // --------------------------------------------------------------------------
  // Test 3: create a Purkinje Network mesh
  // --------------------------------------------------------------------------
 
  cout << endl;
  cout << "TEST 3" << endl;
  cout << "Creating a Purkinje network mesh" << endl;

  /*
  PurkinjeMesh pm;
  pm.read("/data/meshes/heart/purkinje/segm3");
  pm.setup();
  cout << pm << endl;
  */

  return 0;
}
