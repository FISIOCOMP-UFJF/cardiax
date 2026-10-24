#ifndef PIECEWISE_CIRCULATION_HPP
#define PIECEWISE_CIRCULATION_HPP

#include <string>
#include "circulation_model.hpp"
#include "piecewise_pvloop.hpp"

//! PV loop POR FASES (Wang et al. 2025, eLife 14:RP106555), -pw 1.
//!
//! A pressao sai de uma expressao explicita (PiecewisePVLoop::advance_pressure)
//! e vai direto para a mecanica, em vez de ser a incognita de um Newton:
//!   - custo: UMA solucao 3D por passo (contra 4-8 no -circ);
//!   - as fases isovolumetricas sao mantidas por penalidade, nao por restricao;
//!   - circulacao ABERTA (Windkessel de 2 elementos por ventriculo): sem retorno
//!     venoso, sem volume total a conservar, sem Frank-Starling entre batimentos.
//! As fases, os parametros e a estabilizacao estao em piecewise_pvloop.hpp.
class PiecewiseCirculation : public CirculationModel
{
public:
  const char * name() const override { return "PV loop por fases (Wang et al. 2025)"; }

  void configure(const std::string & mshfile, const TimeBase & tb) override;
  void announce(const StepContext & s) override;
  void advance(const StepContext & s) override;
  void end() override;

protected:
  void do_begin() override;

private:
  std::string     base_;
  PiecewisePVLoop pw_;
  int             solves_        = 0;
  int             beat_reported_ = -1;
};

#endif