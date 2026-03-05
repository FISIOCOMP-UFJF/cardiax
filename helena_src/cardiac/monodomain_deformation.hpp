#ifndef MONODOMAIN_DEFORMATION_H
#define MONODOMAIN_DEFORMATION_H

#include "cardiacproblem.hpp"
#include "monodomain.hpp"
#include "../pdes/mech_utils.hpp"

/** 
    This class implements the Monodomain model using a Lagrangian formulation
    to take into account the effects of deformation.    
*/

class MonodomainDeformation : public Monodomain
{
public:

  //! Constitutive assumption for conductivity tensor
  enum CondTensorType { S_ISOTROPIC,    // spatial isotropy
			S_TRANSVERSE,   // spatial transverse isotropy
			S_ORTHOTROPIC,  // spatial orthotropy
			M_ISOTROPIC,    // material isotropy
			M_TRANSVERSE,   // material transverse isotropy
			M_ORTHOTROPIC}; // material orthotropy
  
  //! Default constructor
  MonodomainDeformation();

  //! Default destructor
  ~MonodomainDeformation() {};

  //! Advance one step in time
  void advance(const ArrayMat33 & vec_F);
  
  //! Change conductivity tensor
  void set_conductivity(int cond);

  void init();

private:

  //! Conductivity type
  CondTensorType condtype;

  //! Assemble matrices using deformation gradient F
  void assemble_matrices(const ArrayMat33 & vec_F);

  //! Calculate element stiffness matrix
  void calc_elmat_stiff(const int eindex, const FiniteElement & fe,
                        arma::mat & elmat, const ArrayMat33 & vecF);

  //! Solve the system of linear equations of the parabolic problem
  void solve_parabolic(const ArrayMat33 & vec_F);

  //! Re-calculate the conductivity
  void update_conductivity(const arma::mat33 & F, const arma::vec3 & f0,
                           const arma::vec3 & s0, const arma::vec3 & n0,
                           arma::mat & sigma);
};

#endif


