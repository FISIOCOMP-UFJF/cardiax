#ifndef LAPLACE_HPP
#define LAPLACE_HPP

#include "poisson.hpp"

/** Subclass of Poisson problem to solve the Laplace problem
    \[
      -\Delta u = 0
    \]
*/

class Laplace : public Poisson
{
public:
  Laplace(){};
  ~Laplace(){};

private:
  virtual void assemble_system();
};

#endif
