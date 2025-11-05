#ifndef FUNCTION_HPP
#define FUNCTION_HPP


/** Abstract base class for scalar functions - implemented as functors */

template <typename T> 
class ScalarFunction
{
public:
  virtual ~ScalarFunction() {}
  virtual T operator()(const arma::vec3 & pt) const = 0;

};

/** Concrete class that implements a basic zero function */

template <typename T>
class ZeroFunction : public ScalarFunction<T>
{
public:
  T operator()(const arma::vec3 & pt) const
  {
    return 0.0;
  }
};



#endif
