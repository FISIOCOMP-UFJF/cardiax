#ifndef LEGACY_CIRCULATION_HPP
#define LEGACY_CIRCULATION_HPP

#include <string>
#include <vector>
#include "circulation_model.hpp"

//! Caminho legado: o VD segue a curva de pressao tabelada no <pvloop> da
//! malha, e a pressao do VE sai de um Windkessel com atrio esquerdo
//! (Adams-Bashforth de 4a ordem + Newton sobre o volume da cavidade).
class LegacyCirculation : public CirculationModel
{
public:
  const char * name() const override
  { return "Circulacao legada (curva do <pvloop> + Windkessel)"; }

  void configure(const std::string & mshfile, const TimeBase & tb) override;
  void announce(const StepContext & s) override;
  void advance(const StepContext & s) override;
  void end() override;

  void save_state(State & st) const override;
  void restore_state(const State & st) override;

protected:
  void do_begin() override;

private:
  std::string base_;

  // parametros do <pvloop>
  double c_art_ = 0, c_ven_ = 0, Rmv_ = 0, Rven_ = 0, Rao_ = 0, Rper_ = 0;
  double V_art_zero_ = 0, V_ven_zero_ = 0;
  double E_es_LA_ = 0, A_LA_ = 0, B_LA_ = 0, Tmax_ = 0, tau_ = 0;
  double part_ = 0, pven_ = 0;

  // curva tabelada e historicos (um elemento por passo)
  std::vector<double> curr_time_;                 //!< [ms]
  std::vector<double> p_lv_, p_rv_, p_art_, p_ven_, p_LA_;
  std::vector<double> volume_, volume_rv_;

  // reservatorios do Windkessel no inicio do passo corrente
  double V_LP0_ = 0, V_art0_ = 0, V_ven0_ = 0, V_LA0_ = 10.0;
  double V_LA_zero_ = 10.0;
};

#endif