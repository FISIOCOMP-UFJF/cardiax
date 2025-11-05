#ifndef INTERP_H_
#define INTERP_H_

#include <cmath>
#include <armadillo>

/** Abstract base class for interpolation based on the
    class given in "Numerical Recipes, The art of scientific computing"
*/
struct BaseInterp
{
  int n, mm, jsav;
  int cor; // Search for correlated values
  int dj;
  const double *xx, *yy;

  BaseInterp(arma::vec & x, const double * y, int m)
    : n(x.size()), mm(m), jsav(0), cor(0), xx(x.memptr()), yy(y)
  {
    dj = std::min(1, (int)pow((double)n,0.25));
  }

  virtual ~BaseInterp(){}

  double interp(double x)
  {
    int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
  }

  /** Given x, return j such that x is (insofar as possible) centered
      in the subrange xx[j...j+mm-1], where xx is the stored pointer.
      This routine uses bisection to search.
  */
  int locate(const double x);

  /** Search with correlated values */
  int hunt(const double x);

  /** Interface for interpolation. Concrete base classes must
      implement this member function for interpolation.
  */
  virtual double rawinterp(int jlo, double x) = 0;

};

/** Class for linear interpolation */
struct LinearInterp : BaseInterp
{
  /** Constructor for linear interpolation*/
  LinearInterp(arma::vec & xv, arma::vec & yv)
    : BaseInterp(xv, yv.memptr(), 2) {}

  /** Linear interpolation */
  double rawinterp(int j, double x)
  {
    if (xx[j]==xx[j+1])
      return yy[j]; // table is defective
    else
      return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
  }
};

/** Class for bilinear interpolation on a square grid. */
struct BilinearInterp
{
  int m,n;
  const arma::mat & y;
  LinearInterp x1terp, x2terp;

  /** Constructor */
  BilinearInterp(arma::vec &x1v, arma::vec &x2v, arma::mat &ym)
    : m(x1v.size()), n(x2v.size()), y(ym),
      x1terp(x1v,x1v), x2terp(x2v,x2v) {}

  /** Bilinear interpolation */
  double interp(double x1p, double x2p)
  {
    int i,j;
    double yy, t, u;
    i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
    j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);

    t = (x1p-x1terp.xx[i])/(x1terp.xx[i+1]-x1terp.xx[i]);
    u = (x2p-x2terp.xx[j])/(x2terp.xx[j+1]-x2terp.xx[j]);
    yy = (1.-t)*(1.-u)*y(i,j) + t*(1.-u)*y(i+1,j) + (1.-t)*u*y(i,j+1) + t*u*y(i+1,j+1);
    return yy;
  }
};

#endif
