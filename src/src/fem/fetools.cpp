#include "fetools.hpp"

double FETools::calc_L2_norm(const arma::sp_mat & M, const arma::vec & u)
{
  return sqrt(dot(u,M*u));
}

void FETools::apply_boundary_values(std::map<int,double> nodes,
																		arma::sp_mat & A, arma::vec & b)
{
  typedef std::map <int, double> MapType;
  MapType data = nodes;
  const double p = 1.e+10;

  MapType::const_iterator end = data.end(); 
  for (MapType::const_iterator it = data.begin(); it != end; ++it)
  {
    int index = it->first;
    double value = it->second;
    A(index,index) += p;
    b(index) += p*value;
   }
}

void FETools::apply_boundary_values(std::map<int,double> nodes, 
																		petsc::Matrix & A, petsc::Vector & b)
{
  typedef std::map <int, double> MapType;
  MapType data = nodes;
  const double p = 1.e+10;

  MapType::const_iterator end = data.end(); 
  for (MapType::const_iterator it = data.begin(); it != end; ++it)
  {
    int index = it->first;
    double value = it->second;
    
    A.add(index,index,p);
    b.set(index,value*p);
  }
}

void FETools::apply_boundary_values(std::map<int,double> nodes,
																		petsc::Matrix & A, petsc::Vector & b,
																		petsc::Vector & xbar, bool update_rhs)
{
  typedef std::map <int, double> MapType;
  MapType data = nodes;
  MapType::const_iterator end = data.end();
  int num_rows = nodes.size();
  int * rows = new int[num_rows];
  int k = 0;

  // zero rows and columns and adjust right hand side if necessary
  for (MapType::const_iterator it = data.begin(); it != end; ++it)
    rows[k++] = it->first;

  if (update_rhs)
    A.zero_rows_cols(num_rows, rows, 1.0, xbar, b);
  else
    A.zero_rows_cols(num_rows, rows, 1.0);

  delete [] rows;
}

void FETools::apply_boundary_values(std::map<int,double> nodes, 
                                    petsc::Matrix & A)
{
  typedef std::map <int, double> MapType;
  MapType data = nodes;
  MapType::const_iterator end = data.end();
  int num_rows = nodes.size();
  int * rows = new int[num_rows];
  int k = 0;

  // zero rows and columns and adjust right hand side if necessary
  for (MapType::const_iterator it = data.begin(); it != end; ++it)
    rows[k++] = it->first;

  A.zero_rows_cols(num_rows, rows, 1.0);

  delete [] rows;
}



void FETools::cubic_hermite_deriv(const FESpace & fespace, const Mesh & msh,
																	arma::vec & deriv)
{
	arma::vec elemL;
	elemL.resize(deriv.size());
	
	for(int e=0; e<msh.get_n_elements(); e++)
	{
		std::vector<int> dnums;
		fespace.get_element_dofs(e,dnums);
		arma::vec3 p0 = msh.get_point(dnums[0]);
		arma::vec3 p1 = msh.get_point(dnums[1]);
		double L = p1(0) - p0(0);
		elemL(dnums[0]) = L;
		elemL(dnums[1]) = L;
	}

	for(uint i=0; i<elemL.size(); i++)
		deriv(i) = deriv(i) * (2.0/elemL(i));
}
