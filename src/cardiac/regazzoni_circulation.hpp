#ifndef REGAZZONI_CIRCULATION_HPP
#define REGAZZONI_CIRCULATION_HPP

#include <string>
#include <armadillo>
#include "circulation_model.hpp"
#include "regazzoni2020.hpp"

//! Circulacao 0D fechada (Regazzoni et al. 2022), -circ 1.
//!
//! O 0D e dono dos VOLUMES e pergunta ao 3D, a cada passo, "que pressoes
//! reproduzem estes volumes?". O cardiax aplica pressao (Neumann) e calcula o
//! volume, entao o mapa e invertido numericamente com um Newton 2x2:
//!     F(p_lv, p_rv) = ( V_LV(p) - V_LV_alvo , V_RV(p) - V_RV_alvo )
//! Os termos fora da diagonal importam (o septo e compartilhado). O Jacobiano
//! e reaproveitado entre passos e so refeito quando a convergencia degrada.
class RegazzoniCirculation : public CirculationModel
{
public:
  const char * name() const override { return "Circulacao 0D fechada 3D-0D (Regazzoni 2022)"; }

  void configure(const std::string & mshfile, const TimeBase & tb) override;
  void announce(const StepContext & s) override;
  void advance(const StepContext & s) override;
  void end() override;

protected:
  void do_begin() override;

private:
  //! Acha as pressoes [Pa] que reproduzem os volumes-alvo. INVARIANTE: nao
  //! avanca o modelo celular (varias solucoes 3D por passo fariam a tensao
  //! ativa correr a frente do relogio da circulacao).
  void pressures_for_volumes(double Vlv_target, double Vrv_target, double t,
                             double & plv, double & prv);

  //! Jacobiano dV/dp 2x2 por diferencas finitas (3 solucoes 3D)
  void update_pv_jacobian(double t);

  std::string   base_;
  Regazzoni2020 circ_;
  double        dt_ms_          = 1.0;    //!< passo da circulacao [ms]
  double        t_ms_           = 0.0;    //!< relogio da circulacao [ms]
  int           steps_per_beat_ = 1;

  double        p_lv_cur_ = 0.0;          //!< pressao corrente [Pa]
  double        p_rv_cur_ = 0.0;
  arma::mat22   Jpv_;                     //!< dV/dp [mL/Pa]
  bool          Jpv_valid_ = false;

  int           solves_          = 0;     //!< solucoes 3D (contador de custo)
  int           last_newton_its_ = 0;
  double        last_dV_         = 0.0;
  int           last_solves_     = 0;
};

#endif