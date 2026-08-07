#ifndef CARDIAX_PURKINJE_MESH_HPP
#define CARDIAX_PURKINJE_MESH_HPP

#include "mesh.hpp"

class PurkinjeMesh : public Mesh
{
public:

  //! Default constructor
  PurkinjeMesh();

  //! Default destructor
  ~PurkinjeMesh(){};

  //! Return the PMJ nodes
  const std::vector<int> & get_pmj_nodes() const;

  //! Compute all the bifurcations and PMJ nodes in a Purkinje mesh
  void setup();

  //! Set the root node (it is not a PMJ node, no coupling to tissue)
  void set_root(int id) { root_node = id; };

  //! Get the number of bifurcations - must be used after setup
  size_t get_num_bifurcations() { return bifurcation_nodes.size(); }

  //! An iterator for the bifurcation nodes container
  typedef std::map<int, std::vector<int> >::iterator BifucartionIterator;
  typedef std::map<int, std::vector<int> >::const_iterator const_BifucartionIterator;
  BifucartionIterator begin() { return bifurcation_nodes.begin(); }
  const_BifucartionIterator begin() const { return bifurcation_nodes.begin(); }
  BifucartionIterator end() { return bifurcation_nodes.end(); }
  const_BifucartionIterator end() const { return bifurcation_nodes.end(); }

private:

  //! Container to store the indexes of the bifurcation nodes
  //! (nodes sharing more than 2 elements)
  std::map<int, std::vector<int> > bifurcation_nodes;

  //! Container to store the Purkinje Muscle Junctions (PMJ) nodes
  std::vector<int> pmj_nodes;

  //! Root node
  int root_node;

};

#endif //CARDIAX_PURKINJE_MESH_HPP
