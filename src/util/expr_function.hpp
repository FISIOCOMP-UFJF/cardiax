#ifndef EXPR_FUNCTION_HPP
#define EXPR_FUNCTION_HPP

#include <string>
#include <memory>
#include "fem/fem.h"   // ScalarFunction<double>, arma::vec3

//! A ScalarFunction<double> defined by a string expression in x, y, z.
//!
//! The exprtk dependency is hidden behind a pimpl (struct Impl, defined in
//! expr_function.cpp) so that this header stays lightweight: files that
//! include it -- e.g. poisson.hpp and everything that includes THAT -- no
//! longer pull in the ~48k-line exprtk.hpp and no longer pay its compile cost.
//! exprtk is compiled exactly once, in expr_function.cpp.
//!
//! Constants pi and e are available in expressions. The expression is
//! compiled once at construction; operator() just sets x,y,z and evaluates.
class ExprScalarFunction : public ScalarFunction<double>
{
public:
  explicit ExprScalarFunction(const std::string & expr);
  ~ExprScalarFunction();

  double operator()(const arma::vec3 & p) const override;

private:
  struct Impl;                    // defined in expr_function.cpp
  std::unique_ptr<Impl> impl_;
};

#endif
