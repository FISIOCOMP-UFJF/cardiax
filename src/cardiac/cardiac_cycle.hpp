#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include <memory>
#include <string>
#include <vector>
#include <fstream>

#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include "eikonal_reactive.hpp"

#include "circulation_model.hpp"
#include "cell_coupling.hpp"

class CardiacElectromechanic : private CavityMechanics
{
public:
  explicit CardiacElectromechanic(const std::string &epmodel);
  ~CardiacElectromechanic();

  void config(const std::string &basename);

  EikonalReactive &ref() { return ephy; }

  //! Roda a simulacao com o modelo de circulacao escolhido em config()
  void solve();

private:
  // ---- CavityMechanics -------
  void   apply_pressures(double t_ms, double p_lv, double p_rv) override;
  double volume_LV() override;
  double volume_RV() override;
  double volume_total() override;
  bool   verbose() const override { return !quiet_solve; }

  void configure_fibers();
  void configure_mef(const std::string &cellmodel);
  void describe_active_tension_pipeline(double tascale);

  // ---- um passo ----------------------------------------------------------
  //! Avanca a EP n_sub sub-passos (com estimulo LAT, se o modelo celular o
  //! pede) e atualiza Ta. Devolve true se algum no foi estimulado.
  bool advance_cells(int n_sub, const RunPlan &plan);

  //! N = quantos passos da EP cabem em um passo de acoplamento. Confere que o
  //! relogio da EP cobre a rodada inteira (senao o modelo celular para de
  //! integrar no meio, e Eikonal::advance() vira um no-op silencioso).
  int  cell_substeps(const RunPlan &plan);

  void write_frame(int frame, const RunPlan &plan);

  // ---- tensao ativa ------------------------------------------------------
  //! Le Ta do modelo celular para `ta` e deriva ta_applied e ka
  void update_active_tension();

  //! Preenche `ka`, a rigidez ativa nodal do submodelo de Land, usada pela
  //! mecanica para estabilizar a realimentacao de velocidade. Limpa `ka`
  //! (o que desliga a estabilizacao) sempre que nao se aplica.
  void update_active_stiffness();

  //! Grava Vm e Ta por no, no mesmo ritmo do campo de deslocamento
  void save_node_fields(int frame);

  // ---- componentes ------------------------------------------------------
  EikonalReactive   ephy;
  UpdatedLagrangian elas;
  std::unique_ptr<CirculationModel> circulation;
  CellProfile       cell;          
  LandCoupling      land;         
  TimerSection      timer;
  std::string       filename;

  // ---- tensao ativa ------------------------------------------------------
  double T_ref = 1.0;              //!< escala de Ta fenomenologico (<pvloop T_ref>)
  arma::vec ta;                    //!< Ta nodal, ja em unidades de pressao (ta_gain)

  //! Ganho ESCALAR de Ta: conversao do modelo celular para Pa vezes -tascale
  //! (contratilidade global, o botao de calibracao). Um numero so, aplicado a
  //! Ta e a Ka.
  //!
  //! O ta_scale do MATERIAL e outra coisa: um CAMPO por elemento fixado pela
  //! geometria (1 = miocardio, 0 = plug valvar, intermediario = cicatriz).
  //! Um escalar nao o substitui, e o campo nao substitui o escalar sem
  //! reeditar o XML a cada calibracao. Os dois se multiplicam, na montagem.
  double ta_gain = 1.0;

  //! Tensao ativa EFETIVA por elemento, ta_nodal * ta_scale do material.
  //! Usada so para log e saida (a montagem recebe 'ta' e aplica o mapa por
  //! conta propria), entao so e calculada quando alguem vai le-la.
  arma::vec ta_applied;

  //! Mapa do multiplicador ta_scale por elemento (constante no tempo).
  arma::vec ta_scale_map;

  arma::vec ka;                    //!< rigidez ativa por no, unidades do modelo

  // ---- saida e log -------------------------------------------------------
  arma::vec vm_node;               //!< potencial de membrana por no [mV]
  int  out_every = 8;              //!< grava a malha a cada N passos (-prc)
  int  log_every = 25;             //!< relatorio detalhado a cada N passos (-logc)
  int  n_frames  = 1;              //!< frames alocados no writer
  bool quiet_solve = false;        //!< suprime o log por solucao mecanica

  // ---- mechano-electric feedback: fibre stretch --------------------
  // After each mechanical solve, lambda_f = sqrt(I4f) is evaluated per
  // element, averaged onto the nodes, and handed to the cell model, where it
  // drives h(lambda) and the length-dependent ca50 of the Land submodel.
  bool   lambda_coupling = false;  //!< -lam 1 enables it
  bool   lambda_rate_on  = false;  //!< -lamrate 1 also feeds d(lambda)/dt
  double lambda_min_clip = 0.7;    //!< guard rails on the value handed over
  double lambda_max_clip = 1.3;

  //! Guard rail on d(lambda)/dt, in the CELL MODEL's own time unit (1/ms for
  //! ToRORd). The Land distortion equation has
  //!     ZETAS_ss = (A/cds) * lambda_rate = 249 * lambda_rate,
  //! so lambda_rate = -0.004 /ms (4 lengths/s, the model's Vmax) already
  //! drives ZETAS to -1, where (ZETAS+1) changes sign and the active tension
  //! goes NEGATIVE. A negative Ta is a compressive fibre stress and inverts
  //! elements. The default 0.003 keeps ZETAS above -0.75 with the stock
  //! parameters. Set to 0 (-lamratemax 0) to disable the clamp entirely and
  //! recover the previous, unguarded behaviour.
  double lambda_rate_max_clip = 0.003;

  //! Number of INITIAL calls to update_fiber_stretch() during which the rate
  //! is forced to zero (-lamratedelay N, default 0 = previous behaviour).
  //! The passive inflation takes lambda from 1 to ~1.1 over the first few
  //! steps; differencing that transient produces d(lambda)/dt of order
  //! 0.01 /ms, several times the threshold at which the active tension
  //! changes sign, before the mechanics has settled into anything physical.
  //! lambda itself (and therefore ca50 and h(lambda)) stays coupled
  //! throughout -- only the distortion terms ZETAS/ZETAW are held off.
  //! -lamstab 1 (default) enables the Regazzoni-Quarteroni stabilisation of
  //! the velocity-dependent active tension. Set to 0 to recover the plain
  //! segregated scheme, which is unstable whenever the active stiffness
  //! exceeds the passive one -- useful only for reproducing older runs.
  bool stabilize_active = true;

  int  lambda_rate_delay = 0;
  int  lambda_step_count = 0;      //!< calls to update_fiber_stretch() so far

  arma::vec lambda_node;           //!< lambda_f per node, after clipping
  arma::vec lambda_raw;            //!< lambda_f per node, BEFORE clipping
  arma::vec lambda_prev;           //!< previous lambda_raw, for the rate
  arma::vec lambda_rate_node;      //!< d(lambda_f)/dt, cell-model time units
  bool   lambda_prev_valid = false;

  //! Evaluate lambda_f from the current mechanical state and push it into
  //! the cell model. dt_solver is the interval since the previous call, in
  //! ms (the EP solver's unit); it is converted to the cell model's own unit
  //! internally and is only used for the rate.
  void update_fiber_stretch(double dt_solver);
};

#endif