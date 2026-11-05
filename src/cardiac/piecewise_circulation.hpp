#ifndef PIECEWISE_CIRCULATION_HPP
#define PIECEWISE_CIRCULATION_HPP

#include <string>
#include <armadillo>
#include "circulation_model.hpp"
#include "piecewise_pvloop.hpp"

//! PV loop POR FASES, -pw 1: o pulse.cycle (CycleController) em versao
//! quase-estatica.
//!
//! No pulse cada cavidade tem uma incognita de pressao resolvida no MESMO
//! Newton da mecanica (CavityControl), e a fase escolhe a equacao dela:
//! p = P, V(u) = V_alvo ou p = A + B V(u). O cardiax so aplica pressao
//! (Neumann) e mede volume, entao a mesma restricao e imposta por um Newton
//! EXTERNO nas pressoes das cavidades que nao tem pressao prescrita:
//!     F_i(p) = V_i(p) - V_alvo,i          (VOLUME)
//!     F_i(p) = V_i(p) - (p_i - A_i)/B_i   (AFFINE, ejecao)
//! com dV/dp 2x2 por diferencas finitas, atualizado por Broyden e refeito
//! quando a convergencia degrada (como no -circ). Na convergencia o estado e o
//! mesmo da solucao monolitica; a diferenca e so o custo (varias solucoes 3D
//! por passo nas fases com volume ou Windkessel).
class PiecewiseCirculation : public CirculationModel
{
public:
  const char * name() const override { return "PV loop por fases (pulse.cycle, quase-estatico)"; }

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

  //! Impoe as duas restricoes; deixa a malha no estado convergido e devolve
  //! pressoes e volumes nele.
  void solve_constrained(double t, const Constraint c[2], double p[2], double V[2]);

  //! Uma solucao 3D; devolve os volumes
  void mech_solve(double t, const double p[2], double V[2]);

  //! Colunas `free` de dV/dp por diferencas finitas em torno de (p, V)
  void rebuild_jacobian(double t, const double p[2], const double V[2], const bool free[2]);

  std::string     base_;
  PiecewisePVLoop pw_;

  double    p_cur_[2]  = { 0.0, 0.0 };   //!< ultima pressao convergida [Pa]
  arma::mat22 Jpv_;                      //!< dV/dp [mL/Pa]
  bool      col_valid_[2] = { false, false };   //!< colunas de dV/dp utilizaveis

  double    tol_V_     = 1e-3;           //!< tolerancia do residuo [mL] (-pwtol)
  int       max_it_    = 20;             //!< (-pwmaxit)

  int       solves_        = 0;
  int       failures_      = 0;          //!< passos sem convergencia
  int       last_its_      = 0;
  double    last_res_      = 0.0;
  int       beat_reported_ = -1;
};

#endif
