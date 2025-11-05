#ifndef BIDOMAIN_DEFORMATION_HPP
#define	BIDOMAIN_DEFORMATION_HPP

#include "cardiacproblem.hpp"
#include "bidomain.hpp"
#include "../pdes/mech_utils.hpp"

class BidomainDeformation : public Bidomain
{
public:

    //! Constitutive assumption for conductivity tensor
    enum CondTensorType {
        S_ISOTROPIC,    // spatial isotropy
        S_TRANSVERSE,   // spatial transverse isotropy
        S_ORTHOTROPIC,  // spatial orthotropy
        M_ISOTROPIC,    // material isotropy
        M_TRANSVERSE,   // material transverse isotropy
        M_ORTHOTROPIC   // material orthotropy
    };

    //! Default constructor
    BidomainDeformation();

    //! Default destructor
    ~BidomainDeformation() {};

    //! Advance one step in time with update info from deformation gradient F
    void advance(const ArrayMat33 & vec_F);

    //! Change conductivity tensor
    void set_conductivity(int cond);

private:

    //! Conductivity type
    CondTensorType condtype;

    //! Assemble matrices using deformation gradient F
    void assemble_matrices(const ArrayMat33 & vec_F);

    //! Compute element stiffness matrix for parabolic and elliptic problems.
    void calc_elmat_stiff_ie(const int eindex, const FiniteElement & fe,
                             arma::mat & elmat_i, arma::mat & elmat_ie,
                             const ArrayMat33 & vecF);

    //! Solve elliptic part of bidomain
    void solve_elliptic(const ArrayMat33 & vec_F);

    //! Solve the system of linear equations of the parabolic problem
    void solve_parabolic(const ArrayMat33 & vec_F);

    //! Re-calculate the conductivity
    void update_conductivity(const arma::mat33 & F, const arma::vec3 & f0,
                             const arma::vec3 & s0, const arma::vec3 & n0,
                             arma::mat & sigma, arma::mat & sigma_ie);

};

#endif	/* BIDOMAIN_DEFORMATION_HPP */

