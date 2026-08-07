#include "fedata.hpp"

// ----------------------------------------------------------------------------
// FEData 
// ----------------------------------------------------------------------------

FEData::FEData(int iel, FiniteElement * pfe, std::vector<arma::vec3> & x)
    : fe(pfe)
{ 
  // TODO arrumar o mapping
  //map  = fe->get_mapping(iel);
  quad = Quadrature::create(fe->order(), fe->type());
  
  init();
}

FEData::~FEData()
{
 // TODO
}

void FEData::init()
{
  int i;
  int n_intp = quad->get_num_ipoints();
  int n_node = fe->get_ndof();
   
  // allocate memory and initialize
  for(i = 0; i < n_intp; i++)
  {
    arma::vec * p_shp  = new arma::vec(n_node);
    arma::mat * p_dshp = new arma::mat(n_node, n_node);
    p_shp->zeros();
    p_dshp->zeros();
    shp.push_back( p_shp );
    dshp.push_back( p_dshp );
    detjac.push_back( 0 );
  } 
  // TODO
}

void FEData::update()
{
 // TODO
}

// ----------------------------------------------------------------------------
// MxFEData -> FEData for Mixed Elements 
// ----------------------------------------------------------------------------

MxFEData::MxFEData(int i, const MxFE * mf, std::vector<arma::vec3> & x)
  : fe(mf), points(x)
{
  q = Quadrature::create(2*mf->get_order_u()-2, mf->get_type() );
  //map = fe->get_mapping_ptr(i,x);
  init();
}

MxFEData::~MxFEData()
{
  for(uint i=0; i<shape.size(); i++)
  {
    delete shape[i];
    delete shape_deriv[i];
    delete gradn[i];
    delete jacobInv[i];
    delete jacobians[i];
  }
  
  delete q;
  //delete map;
}

void MxFEData::init()
{
  int nQuadPoints = q->get_num_ipoints();
  int nShapeFuncs = fe->get_ndofs()/fe->get_ndim();
  int nDof = nShapeFuncs;
  int nDim = fe->get_ndim();
  arma::vec3 qpt;
  
  // computing shape functions
  for(int i=0; i < nQuadPoints; i++)
  {
    qpt = q->get_point(i);
    arma::vec * shp  = new arma::vec();
    arma::mat * dshp = new arma::mat();
    fe->calc_shape_u (qpt, *shp);
    fe->calc_deriv_shape_u (qpt, *dshp);
    shape.push_back(shp);
    shape_deriv.push_back(dshp);
  }

  // creating mapping within this class
  for(int i=0; i < nQuadPoints; i++)
  {
    jacobians.push_back(new arma::mat(nDim,nDim));
  }

  // computing mapping (jacobian, jac inv, jac det)
  for(int i=0; i < nQuadPoints; i++)
  {
    arma::mat * J = jacobians[i];
    arma::mat * dshp = shape_deriv[i];

    J->zeros();
    for(int i=0; i<nDim; i++){
      for(int j=0; j<nDim; j++){
	for(int k=0; k<nDof; k++){
	  J->at(i,j) += dshp->at(k,j) * points[k][i];
	}
      }
    }
    
    arma::mat * Ji = new arma::mat(arma::inv(*J));    
    jacobInv.push_back( Ji );
    jacobDet.push_back( arma::det(*J) );

    arma::mat * gr = new arma::mat((*dshp)*(*Ji));
    gradn.push_back( gr );
  }
      
}

/*
inline const arma::vec & MxFEData::get_shape(int iq)
{
  return *(shape[iq]);
}

inline const arma::mat & MxFEData::get_deriv_shape(int iq)
{
  return *(shape_deriv[iq]);
}

const arma::mat & MxFEData::get_gradn(int iq)
{
  return *(gradn[iq]);
}

const Real & MxFEData::get_jacobian_det(int iq)
{
  return jacobDet[iq];
}

inline const Real & MxFEData::get_quad_weight(int iq)
{
  return q->get_weight(iq);
}
*/
