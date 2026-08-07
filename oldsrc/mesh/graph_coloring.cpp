#include <vector>
#include <iterator>
#include "graph_coloring.hpp"

void node_in_element(const int node, const arma::umat & elem,
                     std::vector<uint> & vec_elem)
{
  vec_elem.clear();

  // element
  for(uint i=0; i<elem.n_rows; i++)
  {
    // nodes of the element
    for(uint j=0; j<elem.n_cols; j++)
    {
      if(elem(i,j) == (uint)node)
        vec_elem.push_back(i);
    }
  }  
}

arma::umat comm_matrix(const Mesh & mesh)
{
  int npoin = mesh.get_n_points();
  int nelem = mesh.get_n_elements();
  int nnode = mesh.get_nen();

  std::vector<uint> used;
  arma::umat elem(nelem,nnode);
  arma::umat comm(nelem,nelem);
  elem.zeros();
  comm.zeros();

  // copy element connectivity to matrix -> elem
  for(int e=0; e<nelem; e++)
  {
    std::vector<int> edofs;
    mesh.get_element_pt_nums(e, edofs);
    for(int j=0; j<nnode; j++)
      elem(e,j) = edofs[j];
  }

  // create comm matrix
  for(int i=0; i<npoin; i++)
  {
    node_in_element(i, elem, used);
    for (uint j=0; j < used.size(); j++)
    {
      for(uint k = j+1 ; k < used.size(); k++)
      {
        uint jj, kk;
        jj = used[j];
        kk = used[k];
        comm(jj, kk) = 1;
        comm(kk, jj) = 1;
      }
    }
  }

  return comm;
}

arma::umat greedy_coloring(const arma::umat & L, const int nthreads)
{
  int nelem = L.n_rows;
  uint ncolors;

  std::vector<uint> block_color;
  arma::uvec neighbors;
  arma::uvec color_count(nelem);
  arma::uvec C(nelem);
  C.zeros();

  cout << "comm matrix\n" << L << endl;

  C(0) = 1;
  ncolors = 1;
  color_count(0) = 1;

  for(int i=1; i<nelem; i++)
  {
    block_color.clear();
    neighbors = arma::find( L.row(i) == 1);

    // se os vizinhos ja tem cores associadas a eles, entao
    // marca todas essas cores como bloqueadas
    for(uint j=0; j<neighbors.size(); j++)
    {
      if(C(neighbors(j)) != 0)
        block_color.push_back( C(neighbors(j)) );
    }

    //
    //
    for(uint j=0; j<ncolors; j++)
    {
      if(color_count(j) != (uint)nthreads)
      {
        std::vector<uint> tmp;
        for(uint k=0; k<block_color.size(); k++)
        {
          if(block_color[k] == j)
            tmp.push_back(k);
        }

        bool is_free = tmp.size() == 0 ? true : false;
        if(is_free)
        {
          C(i) = j;
          color_count(j) = color_count(j) + 1;
          break;
        }
      }
    }

    //
    //
    if(C(i) == 0)
    {
      ncolors = ncolors + 1;
      C(i) = ncolors;
      color_count(ncolors-1) = 1;    
    }

  }
  
  cout << "OPA FINAL" << endl;
  return C; 

/*
%% Function that colors the elements using the Greedy Algorithm
function [C,NumberOfColors]=Color(L,threads)
C=zeros(size(L,1),1);
C(1)=1;
NumberOfColors=1;
ColorCount(1)=1;
for i=2:size(L,1)
    BlockedColors=[];
    Neighbors=find(L(i,:)==1);
    for j=1:length(Neighbors)
        if C(Neighbors(j))~=0
            BlockedColors=[BlockedColors C(Neighbors(j))];
        end
    end
    for j=1:NumberOfColors
        if ColorCount(j)~=threads
            IsFree=isempty(find(BlockedColors==j));
            if IsFree
                C(i)=j;
                ColorCount(j)=ColorCount(j)+1;
                break;
            end
        end
    end
    if C(i)==0
        NumberOfColors=NumberOfColors+1;
        C(i)=NumberOfColors;
        ColorCount(NumberOfColors)=1;
    end
end

return
  */
}

/*
function [C,L,NumberOfColors]=GreedyColoringGPU(ELEM,NODE,threads,verbose)
% If the function receives 4 arguments, then verbose mode is enabled and
% the communication matrix is printed on screen.

L=CommMatrix(ELEM);
[C,NumberOfColors]=Color(L,threads);
PlotColoring(ELEM,NODE,C,NumberOfColors,threads)
if nargin==4
    DispCommMatrix(L);
end
return
*/
