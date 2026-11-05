#ifndef PIECEWISE_CIRCULATION_HPP
#define PIECEWISE_CIRCULATION_HPP

#include <string>
#include <armadillo>
#include "circulation_model.hpp"
#include "piecewise_pvloop.hpp"

//! Phase-based PV loop (-pw 1), quasi-static, driven by PiecewisePVLoop.
//!
//! The mechanics only takes cavity pressures (Neumann) and returns volumes, so
//! each phase constraint is imposed by an outer Newton on the pressures of
//! the cavities whose pressure is not prescribed:
//!     F_i(p) = V_i(p) - V_target,i          (VOLUME)
//!     F_i(p) = V_i(p) - (p_i - A_i)/B_i     (AFFINE, ejection)
//! The 2x2 dV/dp is built by finite differences, updated by Broyden and
//! rebuilt when convergence degrades.
class PiecewiseCirculation : public CirculationModel
{
public:
  const char * name() const override { return "PV loop por fases (quase-estatico)"; }

  void configure(const std::string & mshfile, const TimeBase & tb) override;
  void announce(const StepContext & s) override;
  void advance(const StepContext & s) override;
  void end() override;

  void save_state(State & st) const override;
  void restore_state(const State & st) override;

protected:
  void do_begin() override;

private:
  typedef PiecewisePVLoop::Constraint Constraint;

  //! Imposes both constraints; leaves the mesh at the converged state and
  //! returns its pressures and volumes.
  void solve_constrained(double t, const Constraint c[2], double p[2], double V[2]);

  //! One 3D solve; returns the volumes
  void mech_solve(double t, const double p[2], double V[2]);

  //! Columns `free` of dV/dp by finite differences around (p, V)
  void rebuild_jacobian(double t, const double p[2], const double V[2], const bool free[2]);

  std::string     base_;
  PiecewisePVLoop pw_;

  double    p_cur_[2]  = { 0.0, 0.0 };   //!< last converged pressure [Pa]
  arma::mat22 Jpv_;                      //!< dV/dp [mL/Pa]
  bool      col_valid_[2] = { false, false };   //!< usable columns of dV/dp

  double    tol_V_     = 1e-3;           //!< residual tolerance [mL] (-pwtol)
  int       max_it_    = 20;             //!< (-pwmaxit)

  int       solves_        = 0;
  int       failures_      = 0;          //!< steps that did not converge
  int       last_its_      = 0;
  double    last_res_      = 0.0;
  int       beat_reported_ = -1;
};

#endif
