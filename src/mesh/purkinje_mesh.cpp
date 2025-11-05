#include "purkinje_mesh.hpp"

PurkinjeMesh::PurkinjeMesh() : root_node(-1)
{
  // do nothing...
}

const std::vector<int> & PurkinjeMesh::get_pmj_nodes() const
{
  return pmj_nodes;
}

void PurkinjeMesh::setup()
{
  std::vector<int> node_count;
  node_count.resize(points.size(),0);

  // count for each node to how many nodes it is connected
  for(int i=0; i<get_n_elements(); i++)
  {
    Element e = get_element(i);
    const std::vector<int> dofs = e.get_pt_nums();
    int x1 = dofs[0];
    int x2 = dofs[1];
    node_count[x1]++;
    node_count[x2]++;
  }

  // find bifurcation nodes
  int pos = 0;
  for(std::vector<int>::iterator it = node_count.begin();
      it != node_count.end(); ++it)
  {
    // if it has more than 2 neighbors, then it's a bifurcation
    if(*it > 2)
    {
      // find neighbors
      std::vector<int> cn;
      for(size_t i=0; i<get_n_elements(); i++)
      {
        Element e = get_element(i);
        const std::vector<int> dofs = e.get_pt_nums();
        int x1 = dofs[0];
        int x2 = dofs[1];

        if(pos == x1) cn.push_back(x2);
        if(pos == x2) cn.push_back(x1);
      }

      std::pair<int, std::vector<int> > b(pos, cn);
      bifurcation_nodes.insert( b );

      // output
      //cout << "bifurcation node = " << pos << endl;
      //cout << "bifurcation neighbors = ";
      //for(int i=0; i<cn.size(); i++)
      //  cout << cn[i] << "\t";
      //cout << endl << endl;
    }

    // if it has exactly 1 neighbor, then it is the root or an PMJ
    if(*it == 1)
    {
      pmj_nodes.push_back(pos);
    }

    // go to the next node
    pos++;
  }

  // remove root node from PMJ nodes
  if(root_node != -1)
  {
    for (std::vector<int>::iterator it = pmj_nodes.begin();
         it != pmj_nodes.end(); ++it) {
      if (*it == root_node) pmj_nodes.erase(it);
    }
  }
  else
  {
    std::cerr << "Error: you need to set root_node before setup() in PurkinjeMesh\n";
    exit(1);
  }

  cout << "Number of PMJ nodes: " << pmj_nodes.size() << endl;
  cout << "Number of bifurcations: " << bifurcation_nodes.size() << endl;

}