#include "petsc_vector.hpp"
#include <iostream>

namespace petsc
{

Vector::Vector() : _vec(0)
{
  // do nothing
}

Vector::Vector(const Vector & a)
{
  Vec aux = a.vec();
  ierr = VecDuplicate(_vec,&aux);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

Vector::Vector(Vec a)
{
  ierr = VecDuplicate(_vec,&a);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

Vector::~Vector()
{
  destroy();
}
 
void Vector::add(const int i, const double s)
{
  ierr = VecSetValue(_vec, i, s, ADD_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::add(const int n, const int * idx, 
		  const double * values)
{
  ierr = VecSetValues(_vec, n, idx, values, ADD_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::add_vec(Vector & x)
{
  int ni = size();
  int *ix = new int[ni];
  double *y;

  if (ni != x.size()){
    std::cerr << "PETSc::Vectors have different sizes.";
    exit(1);
  }

  ierr = VecGetArray(x.vec(), &y);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  for (int i=0; i<ni; i++) ix[i] = i;

  ierr = VecSetValues(_vec, ni, ix, y, ADD_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecRestoreArray(x.vec(), &y);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  delete [] ix;
}

void Vector::assemble()
{
  ierr = VecAssemblyBegin(_vec);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecAssemblyEnd(_vec);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::copy_values(const int n, Vector & x)
{
  // Copy the contents of x, of size ni into the
  // first ni positions of *this

  if(x.size() < n)
  {
    std::cerr << "PETSc::Vector.copy_vec(): source vector is smaller than n.";
    exit(1);
  }

  int N = size();

  if (N < n)
  {
    std::cerr << "PETSc::Vector.copy_vec(): destination vector is smaller.";
    exit(1);
  }

  int *ix = new int[n];
  double *y;

  ierr = VecGetArray(x.vec(), &y);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  for (int i=0; i<n; i++) ix[i] = i;

  ierr = VecSetValues(_vec, n, ix, y, INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecRestoreArray(x.vec(), &y);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  delete [] ix;
}

void Vector::create(int size)
{
  //
  // Only *Sequential* for now
  //
  ierr = VecCreate(PETSC_COMM_WORLD,&_vec);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecSetSizes(_vec,PETSC_DECIDE,size);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecSetFromOptions(_vec);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecSetOption(_vec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecSet(_vec,0);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

// TODO: need to implement a create function that gets a pointer to some region of memory
//       to be used together with Armadillo v.memptr()

void Vector::destroy ()
{
  if (_vec != NULL){
    ierr = VecDestroy(&_vec);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
}

double Vector::dot(Vector & x)
{
  double s;
  ierr = VecDot(_vec, x.vec(), &s);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  return s;
}

double * Vector::get_array()
{
  double * tmp;
  ierr = VecGetArray(_vec, &tmp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  return tmp;
}

void Vector::save_array(char filename[])
{
  double * tmp;
  ierr = VecGetArray(_vec, &tmp);

  std::ofstream file;
  file.open(filename);

  for(int i=0; i<size(); i++) file << std::scientific << tmp[i] << "\n";

  file.close();

  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::get_data(double *u)
{
  double * tmp;

  ierr = VecGetArray(_vec, &tmp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  for(int i=0; i<size(); i++) u[i] = tmp[i];

  ierr = VecRestoreArray(_vec, &tmp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::get_data(double *u, int k1, int k2)
{
  double * tmp;

  ierr = VecGetArray(_vec, &tmp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  for(int i=0; i<(k2-k1); i++){
    int j = k2-k1+i;
    u[i] = tmp[j];
  }

  ierr = VecRestoreArray(_vec, &tmp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

double Vector::l2norm()
{
  double r;
  ierr = VecNorm(_vec, NORM_2, &r);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  return r;
}

Vector & Vector::operator=(const double s)
{
  assert(s==0);

  ierr = VecZeroEntries (_vec);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  return *this;
}

void Vector::place_array(double * values)
{
  assert(values != NULL);
  ierr = VecPlaceArray(_vec, values);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::reset_array()
{
  ierr = VecResetArray(_vec);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}
 
void Vector::restore_array(double * aa)
{
  ierr = VecRestoreArray(_vec, &aa);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::set(const int i, double s)
{
  ierr = VecSetValue(_vec, i, s, INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::set(const int n, const int * idx, const double * values)
{
  ierr = VecSetValues(_vec, n, idx, values, INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::set(const double s)
{
  ierr = VecSet(_vec, s);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::set_data(const double *u)
{
  double * tmp;

  ierr = VecGetArray(_vec, &tmp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  for(int i=0; i<size(); i++)
    tmp[i] = u[i];

  ierr = VecRestoreArray(_vec, &tmp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::scale (double alpha)
{
  ierr = VecScale(_vec, alpha);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Vector::save_binary(int i)
{
  char filename[20];
  sprintf(filename,"rhsvec_%06d.dat",i);

  PetscViewer viewer;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = VecView(_vec,viewer);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = PetscViewerDestroy(&viewer);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

int Vector::size() const
{
  int s;
  VecGetSize(_vec,&s);
  return s;
}

void Vector::view()
{
  std::cout << std::scientific;
  VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
}

}
