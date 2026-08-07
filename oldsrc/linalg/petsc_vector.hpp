#ifndef PETSC_VECTOR_HPP
#define PETSC_VECTOR_HPP

#include <iostream>
#include <fstream>
#include <cassert>
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

namespace petsc
{

class Vector
{
public:

  Vector();
  Vector(Vec a);
  Vector(const Vector &v);
  ~Vector();
  
  void add(const int i, const double value);
  void add(const int n, const int * idx, const double * values);
  void add_vec(Vector & x);
  void assemble();
  void copy_values(const int n, Vector & x);
  void create(int size);
  void destroy();
  double dot(Vector & y);
  double * get_array();
  void save_array(char filename[]);
  void get_data(double *u);
  void get_data(double *u, int k1, int k2);  
  double l2norm();
  Vector & operator=(const double value);
  void place_array(double * values);
  void reset_array();
  void restore_array(double *);
  void save_binary(int i);
  int size() const;
  void scale(double alpha);
  void set(const int i, double s);
  void set(const int n, const int * idx, const double * values);
  void set(const double s);
  void set_data(const double *u);
  Vec vec() const { return _vec; }
  void view();

private:

  Vec _vec;
  PetscErrorCode ierr;

};


} // namespace PETSc

#endif
