#include "quadrature.hpp"

Quadrature * Quadrature::create(int order, ElementType etype)
{
  if (etype == ELEM_SEGM)
  {
    Quadrature1d * q1d = new Quadrature1d(order);
    return q1d;
  }
  else if (etype == ELEM_TRIG)
  {
    QuadratureTri * qtr = new QuadratureTri(order);
    return qtr;
  }
  else if (etype == ELEM_QUAD)
  {
    QuadratureQuad * qqd = new QuadratureQuad(order);
    return qqd;
  }
  else if (etype == ELEM_TETRA)
  {
    QuadratureTetra * qqd = new QuadratureTetra(order);
    return qqd;
  }
  else if (etype == ELEM_HEXA)
  {
    QuadratureHexa * qqd = new QuadratureHexa(order);
    return qqd;
  }
  else
    throw std::runtime_error("quadrature rule not found");

  return NULL;
}

// -----------------------------------------------------------------------------

Quadrature1d::Quadrature1d(int order) : Quadrature(order)
{
  // 1D quadrature points
  static double quad1d_o1_pt0[3] = {0.0, 0.0, 0.0};
  
  static double quad1d_o2_pt0[3] = {-1.0/sqrt(3), 0.0, 0.0};
  static double quad1d_o2_pt1[3] = {1.0/sqrt(3), 0.0, 0.0};
  
  static double quad1d_n3_pt0[3] = {0.0, 0.0, 0.0};
  static double quad1d_n3_pt1[3] = {-sqrt(3.0/5.0), 0.0, 0.0};
  static double quad1d_n3_pt2[3] = { sqrt(3.0/5.0), 0.0, 0.0};

  switch(order)
  {
  case 0:
  case 1:
    weights.push_back(2.0);
    ipoints.resize(1);
    ipoints[0] = arma::vec3(quad1d_o1_pt0);
    break;
  case 2:
  case 3:
    weights.push_back(1.0);
    weights.push_back(1.0);
    ipoints.resize(2);
    ipoints[0] = arma::vec3(quad1d_o2_pt0);
    ipoints[1] = arma::vec3(quad1d_o2_pt1);
    break;
  case 4:
  case 5:
  case 6:
  case 7:
    weights.push_back(8.0/9.0);
    weights.push_back(5.0/9.0);
    weights.push_back(5.0/9.0);
    ipoints.resize(3);
    ipoints[0] = arma::vec3(quad1d_n3_pt0);
    ipoints[1] = arma::vec3(quad1d_n3_pt1);
    ipoints[2] = arma::vec3(quad1d_n3_pt2);
    break;
  default:
    throw std::runtime_error("invalid 1d integration order");
    break;
  }   
}

// -----------------------------------------------------------------------------

QuadratureTri::QuadratureTri(int order) : Quadrature(order)
{
  // 2D quadrature points - tri
  static double quadtri_o0_pt0[3] = {1./3., 1./3., 0.0};
  
  static double quadtri_o2_pt0[3] = {2./3., 1./6., 0.0};
  static double quadtri_o2_pt1[3] = {1./6., 2./3., 0.0};
  static double quadtri_o2_pt2[3] = {1./6., 1./6., 0.0};
  
  switch(order)
  {
  case 0:
  case 1:
    // weights
    weights.push_back(1.0/2.0);
    // points
    ipoints.resize(1);
    ipoints[0] = arma::vec3(quadtri_o0_pt0);
    break;

  case 2:
    // weights
    for(int i=0; i<3; i++)
      weights.push_back(1.0/6.0);
    // points
    ipoints.resize(3);
    ipoints[0] = arma::vec3(quadtri_o2_pt0);
    ipoints[1] = arma::vec3(quadtri_o2_pt1);
    ipoints[2] = arma::vec3(quadtri_o2_pt2);
    break;
  default:
    throw std::runtime_error("invalid tri integration order");
    break;
  }
}

// -----------------------------------------------------------------------------

QuadratureQuad::QuadratureQuad(int order) : Quadrature(order)
{
  // 2D quadrature data - quad
  const double b = 1.0/sqrt(3);  
  static double quadqd_o0_pt0[3] = {-b, -b, 0.0};
  static double quadqd_o0_pt1[3] = { b, -b, 0.0};
  static double quadqd_o0_pt2[3] = { b,  b, 0.0};
  static double quadqd_o0_pt3[3] = {-b,  b, 0.0};

  const double c = 0.7745966692414834;  
  static double q2p1[3] = {-c,-c, 0};
  static double q2p2[3] = {-c, 0, 0};
  static double q2p3[3] = {-c, c, 0};
  static double q2p4[3] = {0, -c, 0};
  static double q2p5[3] = {0,  0, 0};
  static double q2p6[3] = {0,  c, 0};
  static double q2p7[3] = {c, -c, 0};
  static double q2p8[3] = {c,  0, 0};
  static double q2p9[3] = {c,  c, 0};
    
  switch(order)
  {
  case 0: case 1: case 2:
    // weights
    for(int i=0; i<4; i++)
      weights.push_back(1.0);
    // points
    ipoints.resize(4);
    ipoints[0] = arma::vec3(quadqd_o0_pt0);
    ipoints[1] = arma::vec3(quadqd_o0_pt1);
    ipoints[2] = arma::vec3(quadqd_o0_pt2);
    ipoints[3] = arma::vec3(quadqd_o0_pt3);
    break;
  case 3:
    // weights
    weights.push_back(0.30864197530864212);
    weights.push_back(0.49382716049382724);
    weights.push_back(0.30864197530864212);
    weights.push_back(0.49382716049382724);
    weights.push_back(0.79012345679012341);
    weights.push_back(0.49382716049382724);
    weights.push_back(0.30864197530864212);
    weights.push_back(0.49382716049382724);
    weights.push_back(0.30864197530864212);
    // points
    ipoints.resize(9);   
    ipoints[0] = arma::vec3(q2p1);
    ipoints[1] = arma::vec3(q2p2);
    ipoints[2] = arma::vec3(q2p3);
    ipoints[3] = arma::vec3(q2p4);
    ipoints[4] = arma::vec3(q2p5);    
    ipoints[5] = arma::vec3(q2p6);
    ipoints[6] = arma::vec3(q2p7);
    ipoints[7] = arma::vec3(q2p8);
    ipoints[8] = arma::vec3(q2p9);
    break;    
  default:
    throw std::runtime_error("invalid quad integration order");
    break;
  }
}

// -----------------------------------------------------------------------------

QuadratureTetra::QuadratureTetra(int order) : Quadrature(order)
{
  // 3D quadrature data - tet
  const double a = 0.25;
  const double b = 0.5854101966249685;
  const double c = 0.1381966011250105;  
  static double quadtet_o0_pt0[3] = {a, a, a};      
  static double quadtet_o2_pt0[3] = {b, c, c};
  static double quadtet_o2_pt1[3] = {c, b, c};
  static double quadtet_o2_pt2[3] = {c, c, b};
  static double quadtet_o2_pt3[3] = {c, c, c};
  
  switch(order)
  {
  case 0: case 1:
  {      
    weights.push_back(1./6.);
    ipoints.resize(1);
    ipoints[0] = arma::vec3(quadtet_o0_pt0);
    break;
  }  
  case 2: case 3: case 4:
  {
    const double c = 1.0/24.0;
    for(int i=0; i<4; i++)
      weights.push_back(c);

    ipoints.resize(4);
    ipoints[0] = arma::vec3(quadtet_o2_pt0);
    ipoints[1] = arma::vec3(quadtet_o2_pt1);
    ipoints[2] = arma::vec3(quadtet_o2_pt2);
    ipoints[3] = arma::vec3(quadtet_o2_pt3);
    break;
  }  
  default:
    throw std::runtime_error("invalid tetra integration order");
    break;
  }
}

// -----------------------------------------------------------------------------

QuadratureHexa::QuadratureHexa(int order) : Quadrature(order)
{
  // 3D quadrature data - hex
  const double c = 1.0/sqrt(3);  
  static double quadhx_o0_pt0[3] = {-c, -c, -c};
  static double quadhx_o0_pt1[3] = { c, -c, -c};
  static double quadhx_o0_pt2[3] = { c,  c, -c};
  static double quadhx_o0_pt3[3] = {-c,  c, -c};
  static double quadhx_o0_pt4[3] = {-c, -c,  c};
  static double quadhx_o0_pt5[3] = { c, -c,  c};
  static double quadhx_o0_pt6[3] = { c,  c,  c};
  static double quadhx_o0_pt7[3] = {-c,  c,  c};
  static double quadhx_center[3] = {0 ,  0,  0};
  
  switch(order)
  {
    case -1:
      weights.push_back(8.0);
      ipoints.resize(1);
      ipoints[0] = arma::vec3(quadhx_center);
      break;  
      
    case 0: 
    case 1:
    case 2:
      // weights
      for(int i=0; i<8; i++)
        weights.push_back(1.0);

      // points    
      ipoints.resize(8);
      ipoints[0] = arma::vec3(quadhx_o0_pt0);
      ipoints[1] = arma::vec3(quadhx_o0_pt1);
      ipoints[2] = arma::vec3(quadhx_o0_pt2);
      ipoints[3] = arma::vec3(quadhx_o0_pt3);
      ipoints[4] = arma::vec3(quadhx_o0_pt4);
      ipoints[5] = arma::vec3(quadhx_o0_pt5);
      ipoints[6] = arma::vec3(quadhx_o0_pt6);
      ipoints[7] = arma::vec3(quadhx_o0_pt7);
      break;      
    default:
      throw std::runtime_error("invalid hex integration order");
      break;
  }
}

