#ifndef PIECEWISE_CIRCULATION_HPP
#define PIECEWISE_CIRCULATION_HPP

#include <string>
#include "circulation_model.hpp"
#include "piecewise_pvloop.hpp"

//! PV loop POR FASES (Wang et al. 2025, eLife 14:RP106555 -- a formulacao do
//! Alya), -pw 1.
//!
//! A diferenca estrutural em relacao ao RegazzoniCirculation esta em UMA
//! linha do passo: aqui a pressao sai de uma expressao explicita
//! (PiecewisePVLoop::advance_pressure) e vai direto para a mecanica, em vez
//! de ser a incognita de um Newton 2x2 sobre o volume.
//!
//!   - custo: UMA solucao 3D por passo, contra 4-8 do caminho -circ;
//!   - estabilidade: o acoplamento p <-> V e explicito, entao as fases
//!     isovolumetricas sao mantidas por penalidade e nao por restricao. Ganho
//!     alto demais nas penalidades oscila; baixo demais deixa o volume
//!     derivar. Ver piecewise_pvloop.hpp;
//!   - circulacao ABERTA: um Windkessel de dois elementos por ventriculo, sem
//!     retorno venoso nem atrios. Nao ha volume sanguineo total a conservar,
//!     mas tambem nao ha Frank-Starling entre batimentos.
//!
//! O passo do controlador (-dtp) e dado em ms, independente de -tunit: o
//! artigo e a Tabela 1 estao todos em ms.
class PiecewiseCirculation : public CirculationModel
{
public:
  const char * name() const override
  { return "PV loop por fases (Wang et al. 2025)"; }

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