#include "gmshio.hpp"

namespace {

// hold information about gmsh element types
struct GmshElement
{
  uint id;
  uint dim;
  uint nodes;
};

// create a map to associate our elements with Gmsh element type
std::map<uint, ElementType> elemap;

// elm-type
//  defines the geometrical type of the n-th element:
//  1 2-node line.
//  2 3-node triangle.
//  3 4-node quadrangle.
//  4 4-node tetrahedron.
//  5 8-node hexahedron.
void init_elemap()
{
  elemap[1] = ELEM_SEGM;
  elemap[2] = ELEM_TRIG;
  elemap[3] = ELEM_QUAD;
  elemap[4] = ELEM_TETRA;
  elemap[5] = ELEM_HEXA;
}

}

void GmshIO::read(const std::string & name)
{
  std::ifstream in (name.c_str());
  read_mesh(in);
}

void GmshIO::read_mesh(std::istream& in)
{
  const int buflen = 256;
  char buf[buflen+1];
  int mesh_dim = 3;
  int format=0, size=0;
  double version = 1.0;
  double yacc = 0.0;
  double zacc = 0.0;
  double z0,z1;
  bool zdif = false;

  uint num_elem;

  init_elemap();

  while (!in.eof())
  {
    in >> buf;

    // read header block
    if (std::strncmp(buf,"$MeshFormat",11) == 0)
    {
      in >> version >> format >> size;
      if ((version != 2.0) && (version != 2.1) && (version != 2.2)) {
        std::cerr << " Error: Wrong gmsh file version " << version << "\n";
        exit(1);
      }
      if(format){
        std::cerr << " Error: Unknown data format for mesh\n";
        exit(1);
      }
      // read $MeshFormatEnd
      in >> buf;
    }

    // read nodes block
    else if (std::strncmp(buf,"$Nodes",6) == 0)
    {
      uint num_nodes, node_id;
      in >> num_nodes;
      mesh.set_npoints(num_nodes);
 
      mesh.reserve_points(num_nodes);
      arma::vec3 pt;
      for(uint i=0; i<num_nodes; i++)
      {
        in >> node_id >> pt[0] >> pt[1] >> pt[2];

        if(i==0)
           z0 = pt[2];

        z1 = pt[2];
        if (z0 != z1) zdif = true;

        yacc += pt[1];
        zacc += pt[2];
        mesh.add_point(pt);
      }

      // read $NodesEnd
      in >> buf;

      // find out mesh dimension 1D or 2D or 3D
      if (zdif)
        mesh_dim = 3;
      else if (zacc == 0.0 && yacc != 0.0)
        mesh_dim = 2;
      else if (yacc == 0.0 && zacc == 0.0)
        mesh_dim = 1;

      mesh.set_ndim(mesh_dim);
    }

    // read element data block
    else if(std::strncmp(buf, "$ElementData", 12) == 0)
    {
      int aux, tmp, dpel, numel;
      double val;
      std::string data_name;

      in >> tmp;
      in >> data_name;

      // do nothing. consumes input
      in >> aux;
      in >> val;

      in >> aux;
      in >> tmp >> dpel >> numel;

      for(int i=0; i<numel; i++)
      {
        in >> tmp >> val;
        mesh.add_data(val);
      }
    }

    // read elements block
    else if(std::strncmp(buf, "$Elements", 9) == 0)
    {
      uint id, type, ntags, physical;
      std::string line;

      in >> num_elem;
      mesh.set_nel(num_elem);

      physical = 0;

      for(uint i=0; i<num_elem; i++)
      {
        // read element info
        in >> id >> type >> ntags;

        // read tags
        for(uint j=0; j<ntags; j++)
        {
          uint tag;
          in >> tag;
          if (j==0) physical = tag;
          //else if (j==1) elementary = tag;
        }

        // read connectivity depending on element type
        if (type == 1)
        {
          // build line
          std::vector<int> pnum(2);
          in >> pnum[0] >> pnum[1];

          for(int l=0; l<2; l++) pnum[l]=pnum[l]-1;
	  
          if (mesh_dim == 1)
          {
            mesh.add_elem( Element(ELEM_SEGM,pnum,physical) );
	          mesh.set_nen(2);
	        }
          else if (mesh_dim-1==1)
            mesh.add_boundary_elem( Element(ELEM_SEGM,pnum,physical) );
          else
            print_error("gmshio.cpp", "reading", "1D element in 3D mesh");
        }
        else if (type == 2)
        {
          // build triangle
          std::vector<int> pnum(3);
          in >> pnum[0] >> pnum[1] >> pnum[2];

          for(int l=0; l<3; l++) pnum[l]=pnum[l]-1;

          if (mesh_dim == 2){
            mesh.set_nen(3);
            mesh.add_elem( Element(ELEM_TRIG,pnum,physical) );
          }
          else if (mesh_dim-1==2)
	        {
            mesh.add_boundary_elem( Element(ELEM_TRIG,pnum,physical) );
          }
        }
        else if (type == 3)
        {
          // build quadrangle
          std::vector<int> pnum(4);
          in >> pnum[0] >> pnum[1] >> pnum[2] >> pnum[3];

          for(int l=0; l<4; l++) pnum[l]=pnum[l]-1;

          if (mesh_dim == 2){
            mesh.set_nen(4);
            mesh.add_elem( Element(ELEM_QUAD,pnum,physical) );
          }
          else if (mesh_dim-1==2)
            mesh.add_boundary_elem( Element(ELEM_QUAD,pnum,physical) );
        }
        else if (type == 4)
        {
          // build tetrahedron
          std::vector<int> pnum(4);
          in >> pnum[0] >> pnum[1] >> pnum[2] >> pnum[3];

          for(int l=0; l<4; l++) pnum[l] = pnum[l]-1;

          mesh.set_nen(4);
          mesh.add_elem( Element(ELEM_TETRA,pnum,physical) );
        }
        else if (type == 5)
        {
          // build hexahedron
          std::vector<int> pnum(8);
          in >> pnum[0] >> pnum[1] >> pnum[2] >> pnum[3];
          in >> pnum[4] >> pnum[5] >> pnum[6] >> pnum[7];

          for(int l=0; l<8; l++) pnum[l] = pnum[l]-1;

          mesh.set_nen(8);
          mesh.add_elem( Element(ELEM_HEXA,pnum,physical) );
        }
        else
        {
	        throw std::runtime_error("gmsh element not supported");
        }
      }

    }
  }

  // if triangular mesh, check for positive volumes
  mesh.set_tri_positive_volume();

  // if tetrahedral mesh, check for positive volumes
  mesh.set_tet_positive_volume();

}
