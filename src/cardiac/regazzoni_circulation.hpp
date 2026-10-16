#ifndef REGAZZONI_CIRCULATION_HPP
#define REGAZZONI_CIRCULATION_HPP

#include <string>
#include <armadillo>
#include "circulation_model.hpp"
#include "regazzoni2020.hpp"

//! Circulacao 0D fechada (Regazzoni et al. 2022), -circ 1.
//!
//! O modelo 0D e dono dos VOLUMES das cavidades e pergunta ao 3D, a cada
//! passo, "que pressoes reproduzem estes volumes?". O fenicsx-pulse responde
//! isso direto (um multiplicador de Lagrange). O cardiax aplica pressao como
//! condicao de Neumann e calcula o volume, entao o mapa corre no sentido
//! oposto e e invertido numericamente:
//!
//!     F(p_lv, p_rv) = ( V_LV(p) - V_LV_alvo , V_RV(p) - V_RV_alvo )
//!
//! com Newton 2x2. Os termos fora da diagonal importam: o septo e
//! compartilhado, entao subir a pressao do VE muda o volume do VD. O
//! Jacobiano e reaproveitado entre passos e so refeito quando a convergencia
//! degrada, o que mantem o custo perto de uma solucao 3D por iteracao.
class RegazzoniCirculation : public CirculationModel
{
public:
  const char * name() const override
  { return "Circulacao 0D fechada 3D-0D (Regazzoni 2022)"; }

  void configure(const std::string & mshfile, const TimeBase & tb) override;
  void announce(const StepContext & s) override;
  void advance(const StepContext & s) override;
  void end() override;

protected:
  void do_begin() override;

private:
  //! Inverte o modelo 3D: acha as pressoes que reproduzem os volumes-alvo.
  //! INVARIANTE: nao avanca o modelo celular. Ele roda varias solucoes 3D por
  //! passo, e avancar a EP aqui faria a tensao ativa correr a frente do
  //! relogio da circulacao por um fator igual ao numero de solucoes.
  void pressures_for_volumes(double Vlv_target, double Vrv_target, double t,
                             double & plv_mmHg, double & prv_mmHg);

  //! Reconstroi o Jacobiano dV/dp 2x2 por diferencas finitas (3 solucoes)
  void update_pv_jacobian(double t);

  std::string   base_;
  Regazzoni2020 circ_;
  double        dt_s_            = 1.0e-3;  //!< passo da circulacao [s]
  double        t_s_             = 0.0;     //!< relogio da circulacao [s]
  int           steps_per_beat_  = 1;
  double        mmHg_to_p_       = 133.322387415;

  double        p_lv_cur_ = 0.0;            //!< pressao corrente do VE [unid. do modelo]
  double        p_rv_cur_ = 0.0;
  arma::mat22   Jpv_;                       //!< dV/dp [mL / unid. de pressao]
  bool          Jpv_valid_ = false;

  int           solves_          = 0;       //!< solucoes 3D feitas (contador de custo)
  int           last_newton_its_ = 0;
  double        last_dV_         = 0.0;
  int           last_solves_     = 0;
};

#endif