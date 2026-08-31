#ifndef FEDATA_HPP
#define FEDATA_HPP

#include "fe.hpp"
#include "mixed_fe.hpp"
#include "mapping.hpp"
#include "quadrature.hpp"

/*!
 *  Facade to integrate the following classes:
 *  - Mapping
 *  - Quadrature
 *  - FiniteElement
*/
class FEData
{
public:
    
  FEData(int iel, FiniteElement * pfe, std::vector<arma::vec3> & x);
  ~FEData();
  
  void init();
  void update();
  
  const double & dx();
  const arma::vec & shape();
  const arma::mat & dshape();  
  
private:
  
  FiniteElement * fe;
  Mapping * map;
  Quadrature * quad;
  
  //! Store FE information on each integration point
  std::vector<arma::vec*> shp;
  std::vector<arma::mat*> dshp;
  std::vector<double> detjac; 

};

/*!
 *  Facade to integrate the following classes:
 *  - Mapping
 *  - Quadrature
 *  - MixedFiniteElement
*/
class MxFEData
{
private:

  const MxFE * fe;  
  Quadrature * q;  

  // Element coordinates
  std::vector<arma::vec3> & points;

  // Mapping data
  std::vector<arma::vec*> shape;
  std::vector<arma::mat*> shape_deriv;
  std::vector<arma::mat*> gradn;
  std::vector<arma::mat*> jacobians;
  std::vector<arma::mat*> jacobInv;
  std::vector<double> jacobDet;

public:

  MxFEData(int iel, const MxFE * fe, std::vector<arma::vec3> & x);
  ~MxFEData();

  void init();
  
  const arma::vec & get_shape(int iq)       { return *(shape[iq]); } 
  const arma::mat & get_gradn(int iq)       { return *(gradn[iq]); }
  const arma::mat & get_deriv_shape(int iq) { return *(shape_deriv[iq]); }
  const double    & get_jacobian_det(int iq){ return jacobDet[iq]; }
  const double    & get_quad_weight(int iq) { return q->get_weight(iq); }

};

#endif
