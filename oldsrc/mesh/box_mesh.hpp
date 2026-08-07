#ifndef CARDIAX_BOX_MESH_H
#define CARDIAX_BOX_MESH_H

#include "mesh.hpp"

class BoxMesh : public Mesh
{
public:

  BoxMesh(double x0, double x1, double y0, double y1, double z0, double z1, int n);
  BoxMesh(double x0, double x1, double y0, double y1, double z0, double z1, int nx, int ny, int nz);

private:

  //! Bounding box [xmin,xmax] x [ymin,ymax] x [zmin,zmax]
  double xmin, xmax, ymin, ymax, zmin, zmax;

  //! Number of elements in x direction
  int nx;

  //! Number of element in y direction
  int ny;

  //! Number of element in z direction
  int nz;

  //! Create regular mesh for the bounding box with given spacing
  void create_mesh();

  //! Check if a given point is on the boundary
  bool point_in_boundary(arma::vec3 & pt);

};

#endif //CARDIAX_BOX_MESH_H
