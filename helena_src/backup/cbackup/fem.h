#ifndef _FEM_H_
#define _FEM_H_

/**
 * @file   fem.h
 * @author Bernardo M. Rocha <bernardo@lncc.br>
 * @date   Wed Mar 10 16:24:14 2010
 * 
 * @brief  Finite element routines header file
 * 
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

#include "datatype.h"
#include "util.hpp"
#include "formats.h"
#include "mesh.hpp"

/**
 * Mesh
 */
/*
struct mesh_type {
  int ndim;
  int ndof;
  int nel;
  int nen;
  int nint;
  int nx;
  int ny;
  int eltype; 
  int **localToGlobal; // localToGlobal(nen,nel)
  real **nodes;        // nodes(ndim,ndof)
  char *nodefile; 
  char *elemfile;
  char *stimfile;
  char *axesfile;
};

typedef struct mesh_type mesh;
*/

/**
 * Build global stiffness and mass matrix
 */
void
fem_get_matrices (const Mesh & msh, int *nz, 
		  csr_matrix **csrK, csr_matrix **csrM);

/**
 * Prepare stiffness matrices for solution
 * Create the LHS
 *  - forward Euler   A = (M - dt*K) 
 *  - backward Euler  A = (M + dt*K)
 *  - crank-Nicolson  Al = (M + dt/2 * K)
 *                    Ar = (M - dt/2 * K)
 *  - TODO: add theta
 *     theta = 1   implicit Euler
 *     theta = 0   forward Euler
 *     theta = 1/2 Crank-Nicolson
 *
 * ***Note*** Crank Nicolson in use
 *
 */
void
fem_setup_matrices (int ndof, Real k,
		    csr_matrix *K, csr_matrix *M, csr_matrix **A);

/**
 * Assemble global stiffness
 */
void
fem_assemble_stiffness (const Mesh & msh, csr_matrix *K, int *nz);

/**
 * Assemble global mass
 */
void
fem_assemble_mass (const Mesh & msh, csr_matrix *M, int *nz);

/**
 * Assemble stimulus vector
 */
void
fem_global_stim (int ndim, int nel, int nen, int nint, int eltype, 
		 int **conec, Real **nodes, Real *f, Real *iapp);

/**
 * Assemble load vector for given RHS function
 *
 * @param[in]  nel     Number of elements
 * @param[in]  nen     Number of nodes per element
 * @param[in]  nint    Number of integration points
 * @param[in]  eltype  Element type
 * @param[in]  conec   Element conectivities array
 * @param[in]  nodes   Nodes coordinates array
 * @param[out] f       Global load vector
 * @param[in]  ptRhs   Pointer to RHS function
 */
void
fem_global_rhs (int ndim, int nel, int nen, int nint, int eltype, 
		int **conec, Real **nodes, Real *f, 
		Real (*ptRhs)(Real, Real), Real t);

/**
 * Read mesh file
 */
void
mesh_read (get_mesh * msh, const char *filename, element_list *el);

/**
 * Read elements from .elem file
 */
void
mesh_read_elements (get_mesh *msh);

/**
 * Read nodes from .pts file
 */
void
mesh_read_nodes (get_mesh *msh);

/**
 * Deallocate mesh struct
 */
void
delete_mesh (get_mesh *msh);

/**
 * Deallocate element list
 */
void
delete_element_list (element_list *el);

#endif
