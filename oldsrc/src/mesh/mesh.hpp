#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <map>
#include "fem/element.hpp"
#include "fem/mapping.hpp"

enum MappingType
{
    MAPPING_AFFINE=1,
    MAPPING_ISOPARAMETRIC=2,
};

enum PropertyType
{
    ISOTROPIC,
    TRANSVERSELY_ISOTROPIC,
    ORTHOTROPIC
};

// Forward declaration

class Mapping;
class Mapping;
class SurfaceMapping;
class FiniteElement;

// -----------------------------------------------------------------------------

/** Mesh representation.

    Derived classes implement some common/structured meshes, like
    unit square, unit cube, cable, and so on.

    It supports several input/output formats.
    Supported formats:
       - CARP format (.pts, .elem, .bdn)
       - Gmsh format (.msh)

    To implement:
       - Abaqus (.inp)

    TODO:
       - improve the implementation, implement Mesh as an abstract base class ?
       and then implement Gmsh, CARP and other derived classes for reading
*/
class Mesh
{
public:

  //! Default constructor
  Mesh() : n_dim(), n_points(), n_elements() {}

  //! Destructor
  virtual ~Mesh() {}

  //! Add point to points
  void add_point(const arma::vec3 & pt) { points.push_back(pt); }

  //! Add an element (elem) to the elems std::vector.
  void add_elem(const Element & elem) { elems.push_back(elem); }

  //! Add element elem to boundary_elems
  void add_boundary_elem(const Element & elem) { boundary_elems.push_back(elem); }

  //! Add data to elems
  void add_data(const double & val) { data.push_back(val); }

  //! Allocates and returns pointer to Mapping
  Mapping * element_mapping(int elnum, std::vector<arma::vec3> & x) const;

  //! Return number of dimensions of the mesh
  int get_n_dim() const { return n_dim; }

  //! Return number of nodes per element
  int get_nen() const { return n_node; }

  //! Return number of elements
  int get_n_elements() const { return elems.size(); }

  //! Return number of elements
  int get_n_boundary_elements() const { return boundary_elems.size(); }

  //! Return number of vertices
  uint get_n_points() const { return points.size(); }

  //! Get the i-th element
  const Element& get_element(int i) const { return elems[i]; }

  //! Get the i-th surface element
  const Element& get_boundary_element(int i) const { return boundary_elems[i]; }

  //! Get the i-th point
  const arma::vec3& get_point(int i) const { return points[i]; }

  //! Get reference to elements vector
  const std::vector<Element> & get_elements() const { return elems; }

  //! Get reference to surface elements vector
  const std::vector<Element> & get_selements() const { return boundary_elems; }

  //! Get reference to points
  const std::vector<arma::vec3> & get_points() const { return points; }

  //! Get reference to constrained nodes map
  const std::map<int,double> & get_constrained_nodes() const { return bnodes_map; }

  //! Get node numbers of the element
  void get_element_pt_nums(int elnum, std::vector<int> & pnums) const;

  //! Get Point numbers of surface element
  void get_boundary_element_pt_nums(int selnum, std::vector<int> & pnums) const;

  //! Compute centroid of element
  arma::vec3 get_element_centroid(int elnum) const;

  //! Reference element mapping object that corresponds to element i.
  Mapping get_element_mapping(int i) const;

  //! Reference element mapping object that corresponds to element i.
  Mapping get_element_mapping(int elnum, std::vector<arma::vec3> & x) const;

  //! Ref boundary elem mapping object that corresponds to boundary elem i
  SurfaceMapping get_boundary_element_mapping(int i) const;

  //! Ref boundary elem mapping object that corresponds
  //! to boundary element i with given points
  SurfaceMapping get_boundary_element_mapping(int i, std::vector<arma::vec3> & x) const;

  //! Get reference to data
  const double & get_data(int i) const { return data[i]; }

  //! Return data vector size
  uint get_data_size() const { return data.size(); }

  //! Get region index of element i
  int get_element_index(int i) const { return elems[i].get_index(); }

  //! Get property type
  PropertyType get_prop_type() const { return proptype; }

  //! Read mesh files given basename
  void read(const std::string & filename);

  //! Read mesh files given basename
  void read_xml(const std::string & filename);

  //! Read the 3 vectors for fiber orientation
  void read_fiber(istream & in, int code,
                  arma::vec3 & f,
                  arma::vec3 & t,
                  arma::vec3 & n);

  //! Reserve space for std::vector of points.
  void reserve_points(const uint num);

  //! Reserve space for std::vector of elements
  void reserve_elements(const uint num);

  //! Reserve space for std::vector of points.
  //void reserve_points(const uint num) { points.reserve(num); }

  //! Set number of vertices
  void set_npoints(uint np) { n_points = np; }

  //! Set number of elements
  void set_nel(uint ne) { n_elements = ne; }

  //! Set number of elements
  void set_nen(uint nnode) { n_node = nnode; }

  //! Set number of dimensions
  void set_ndim(uint dim) { n_dim = dim; }

  void set_element_fiber (int i, arma::vec3 & f) { elems[i].set_fiber(f); }
  void set_element_trans (int i, arma::vec3 & s) { elems[i].set_trans(s); }
  void set_element_normal(int i, arma::vec3 & n) { elems[i].set_normal(n); }

  //! If triangular mesh, check for positive volumes
  void set_tri_positive_volume();

  //! If tetrahedral mesh, check for positive volumes
  void set_tet_positive_volume();

  //! Change element node ordering
  void swap_nodes(int e, int a, int b);

  //! Change the coordinates of the point i.
  void update_point(const int i, const arma::vec3 & pt);

  //! Friend function to print mesh object information such as
  //! number of nodes, number of elements and so on.
  friend std::ostream & operator<<(std::ostream&, const Mesh &);


  //
  // ITERATORS - Element
  //
  typedef std::vector<Element>::iterator ElementIterator;
  typedef std::vector<Element>::const_iterator const_ElementIterator;

  ElementIterator begin() { return elems.begin(); }
  const_ElementIterator begin() const { return elems.begin(); }

  ElementIterator end() { return elems.end(); }
  const_ElementIterator end() const { return elems.end(); }

protected:

  //! Number of dimensions
  int n_dim;

  //! Number of nodes per element
  int n_node;

  //! Number of nodes
  int n_points;

  //! Number of elements
  int n_elements;

  //! Number of surface elements
  int n_boundary_elements;

  //! Nodes
  std::vector<arma::vec3> points;

  //! List of elements
  std::vector<Element> elems;

  //! List of boundary elements (in 2D lines - in 3D quad/tri)
  std::vector<Element> boundary_elems;

  //! To store data from elements (used in gmsh mesh)
  std::vector<double> data;

  //! Constrained nodes
  std::map<int, double> bnodes_map;

  //! Type of element fiber
  PropertyType proptype;

};

#endif