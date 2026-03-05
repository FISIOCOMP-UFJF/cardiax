#include <stdexcept>

#include "petsc_matrix.hpp"

namespace petsc
{

void Matrix::create(int rows, int cols, int nz)
{
  //
  // Only *Sequential* for now
  //
  //ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, rows, cols, nz, NULL, &_mat);
  //CHKERRABORT(PETSC_COMM_WORLD,ierr);

  //
  // MatAIJ For parallel with SuperLU or MUMPs
  // 
  //ierr = MatCreate(PETSC_COMM_WORLD, &_mat);
  //CHKERRABORT(PETSC_COMM_WORLD,ierr);

  //ierr = MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
  //CHKERRABORT(PETSC_COMM_WORLD,ierr);

  //ierr = MatSetType(_mat, MATMPIAIJ);
  //CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &_mat);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSetSizes(_mat, rows, cols, PETSC_DECIDE, PETSC_DECIDE);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSetType(_mat, MATSEQAIJ);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
//nodal block size
  //ierr = MatSetBlockSize(_mat,3);
  //CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSetUp(_mat);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  
  //ierr = MatMPIAIJSetPreallocation(_mat, nz, NULL, nz, NULL);
  //CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSeqAIJSetPreallocation(_mat, nz, NULL);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr =MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSetFromOptions(_mat);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Matrix::add(const int i, const int j, double s)
{
  ierr = MatSetValue(_mat, i, j, s, ADD_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Matrix::add(const int m, const int n, const int * idxm,
                 const int * idxn, const double * values)
{
  ierr = MatSetValues(_mat, m, idxm, n, idxn, values, ADD_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Matrix::set(const int i, const int j, double s)
{
  ierr = MatSetValue(_mat, i, j, s, INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Matrix::set(const int m, const int n, const int * idxm,
                 const int * idxn, const double * values)
{
  ierr = MatSetValues(_mat, m, idxm, n, idxn, values, INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

Matrix &  Matrix::operator=(const double s)
{
  assert(s==0);

  ierr = MatZeroEntries (_mat);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  return *this;
}

void Matrix::mult (const Vector & x, Vector & y)
{
  ierr = MatMult(_mat, x.vec(), y.vec());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Matrix::assemble ()
{
  ierr = MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

uint Matrix::get_nnz()
{
  uint nnz;
  MatInfo info;
  double nz_used;

  ierr = MatGetInfo(_mat, MAT_LOCAL, &info);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  nz_used = info.nz_used;
  nnz = (uint) nz_used;
  return nnz;
}

void Matrix::get_CSR(int *n, int *ia, int *ja, double *v)
{
  //
  // WARNING-> RETURN 1 based arrays! For usage with PARDISO
  //
  int ln, lnnz;
  PetscBool done;
  const int *pia, *pja;
  double * vals;
  lnnz = get_nnz();

  ierr = MatGetRowIJ(_mat, 1, PETSC_FALSE, PETSC_FALSE, &ln, &pia, &pja, &done);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSeqAIJGetArray(_mat, &vals);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  //std::cout<<"Não zeros: "<<pia[ln]<<"\n";

  // copy data
  *n = ln;
  for(int i=0; i<ln+1; i++) ia[i] = pia[i];
  for(int i=0; i<lnnz; i++) ja[i] = pja[i];
  for(int i=0; i<lnnz; i++) v[i] = vals[i];

  ierr = MatRestoreRowIJ(_mat, 1, PETSC_FALSE, PETSC_FALSE, &ln, &pia, &pja, &done);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSeqAIJRestoreArray(_mat, &vals);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  if(!done) throw std::runtime_error("error in MatGetRowIJ");
}

void Matrix::store_CSR()
{
  //
  // WARNING-> RETURN 1 based arrays! For usage with PARDISO
  //
  int ln, lnnz;
  PetscBool done;
  const int *pia, *pja;
  double * vals;
  lnnz = get_nnz();

  ierr = MatGetRowIJ(_mat, 1, PETSC_FALSE, PETSC_FALSE, &ln, &pia, &pja, &done);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSeqAIJGetArray(_mat, &vals);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  //std::cout<<"Não zeros: "<<pia[ln]<<"\n";

  std::ofstream file;
  //string aux = basename + string("_matrix.txt");
  file.open("csr_matrix.txt");

  // save data
  file << "ia: " << ln+1 << "\n";
  for(int i=0; i<ln+1; i++) file << pia[i] << " ";
  file << "\n";
  file << "ja: " << lnnz << "\n";
  for(int i=0; i<lnnz; i++) file << pja[i] << " ";
  file << "\n";
  file << "vals: " << lnnz << "\n";
  for(int i=0; i<lnnz; i++) file << std::scientific << vals[i] << " ";
  file << "\n";
  file.close();

  ierr = MatRestoreRowIJ(_mat, 1, PETSC_FALSE, PETSC_FALSE, &ln, &pia, &pja, &done);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSeqAIJRestoreArray(_mat, &vals);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  if(!done) throw std::runtime_error("error in MatGetRowIJ");
}

void Matrix::view()
{
  MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
}

void Matrix::set_symmetric()
{
  ierr = MatSetOption(_mat, MAT_SYMMETRIC, PETSC_TRUE);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = MatSetOption(_mat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

uint Matrix::size()
{
  PetscInt nr, nc;

  ierr = MatGetSize(_mat, &nr, &nc);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  if(nr==nc) return nr;
  else return 0;
}

void Matrix::zero_rows_cols(int num_rows, int * rows, double diag)
{
  ierr = MatZeroRowsColumns(_mat, num_rows, rows, 1.0, PETSC_NULL, PETSC_NULL);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Matrix::zero_rows_cols(int num_rows, int * rows, double diag,
                            petsc::Vector & x, petsc::Vector & b)
{
  ierr = MatZeroRowsColumns(_mat, num_rows, rows, 1.0, x.vec(), b.vec());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void Matrix::setNullSpace(petsc::Vector *coord)
{
   MatNullSpace matnull;
   Vec vec_coord;
   vec_coord = coord->vec();
   ierr = VecSetBlockSize(vec_coord,3);

 CHKERRABORT(PETSC_COMM_WORLD,ierr);
 ierr=MatNullSpaceCreateRigidBody(vec_coord,&matnull);
 CHKERRABORT(PETSC_COMM_WORLD,ierr);
 ierr=MatSetNearNullSpace(_mat,matnull);
 CHKERRABORT(PETSC_COMM_WORLD,ierr);
 ierr = MatNullSpaceDestroy(&matnull);
 CHKERRABORT(PETSC_COMM_WORLD,ierr);
  //std::cout<<"\n\nEntrou\n\n"<<std::endl;
}

}



