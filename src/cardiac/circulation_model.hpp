#ifndef CIRCULATION_MODEL_HPP
#define CIRCULATION_MODEL_HPP

#include <fstream>
#include <string>

// ======================================================================
//  Interface entre o laco do cardiac_cycle e os modelos de circulacao
// ======================================================================
//
//  O cardiac_cycle e dono do que e comum aos tres esquemas de acoplamento
//  (celulas, tensao ativa, realimentacao mecano-eletrica, saida). Cada
//  CirculationModel e dono do que muda: de onde vem a pressao das cavidades
//  a cada passo.
//
//    LegacyCirculation     curva de pressao tabelada no <pvloop> + Windkessel
//    RegazzoniCirculation  circulacao 0D fechada (Regazzoni 2022), -circ 1
//    PiecewiseCirculation  PV loop por fases (Wang et al. 2025),   -pw 1
//
//  Todos os tempos desta interface estao em MILISSEGUNDOS. A conversao para
//  a unidade do solver da EP (-tunit) acontece em um unico lugar, no laco.

//! O lado 3D, como o modelo de circulacao o enxerga: aplica pressoes nas
//! cavidades e devolve volumes. E a UNICA porta entre os dois mundos.
class CavityMechanics
{
public:
  virtual ~CavityMechanics() {}

  //! UMA solucao mecanica com as pressoes dadas [unidades de pressao do modelo]
  virtual void apply_pressures(double t_ms, double p_lv, double p_rv) = 0;

  virtual double volume_LV() = 0;            //!< [mL]
  virtual double volume_RV() = 0;            //!< [mL]
  virtual double volume_total() = 0;         //!< cavidade total [mL]

  virtual double pressure_unit_per_mmHg() const = 0;
  virtual double pressure_unit_per_kPa()  const = 0;

  //! false nos passos em que o log detalhado esta suprimido
  virtual bool verbose() const = 0;
};

//! Entradas de tempo ja resolvidas por config(), entregues a configure()
struct TimeBase
{
  double      solver_ms = 1000.0;  //!< ms por unidade de tempo da EP (-tunit)
  std::string unit      = "s";     //!< "s" ou "ms"
  double      dt        = 0.0;     //!< passo da EP [unidades do solver]
  double      ceiling   = -1.0;    //!< -t [unidades do solver]; <= 0 se ausente
  int         beats     = 1;       //!< -beats
};

//! O que o laco precisa saber da rodada. Preenchido por configure().
struct RunPlan
{
  int    beats             = 1;
  double run_ms            = 0.0;   //!< duracao total da rodada
  double step_ms           = 0.0;   //!< passo de acoplamento
  int    nsteps            = 0;     //!< numero de passos
  double period_ms         = 0.0;   //!< periodo cardiaco (0 = nao se aplica)
  double stim_offset_ms    = 0.0;   //!< origem do relogio do estimulo LAT
  int    span_margin_steps = 2;     //!< folga do span da EP, em passos da EP
  int    default_out_every = 8;     //!< default de -prc
  bool   store_stress      = false; //!< gravar tambem o campo de tensoes
};

//! Contexto de um passo, entregue ao modelo
struct StepContext
{
  int    index;    //!< 0-based
  int    nsteps;
  double t_ms;     //!< instante ao FIM do passo
  double dt_ms;
  bool   log;      //!< este passo imprime o relatorio detalhado
};

class CirculationModel
{
public:
  virtual ~CirculationModel() {}

  virtual const char * name() const = 0;

  //! Le as entradas (XML da malha + linha de comando) e preenche o plano.
  virtual void configure(const std::string & mshfile, const TimeBase & tb) = 0;

  const RunPlan & plan() const { return plan_; }

  //! Estado inicial. Chamado depois de pre_solve() e das condicoes iniciais
  //! celulares, com a tensao ativa de t = 0 ja disponivel.
  void begin(CavityMechanics & mech) { mech_ = &mech; do_begin(); }

  //! Chamado no inicio de TODO passo (banners, cabecalho do log)
  virtual void announce(const StepContext &) {}

  //! Um passo: calcula as pressoes das cavidades e resolve a mecanica.
  //! As celulas e a tensao ativa ja foram avancadas pelo laco.
  virtual void advance(const StepContext & s) = 0;

  //! Fim da rodada: resumo e fechamento de arquivos
  virtual void end() {}

protected:
  virtual void do_begin() = 0;

  //! Reduz `beats` ate a rodada caber no teto -t. -t e um TETO, nunca um
  //! alvo: se nao cobrir os batimentos pedidos, reduz o numero deles em vez
  //! de deixar o relogio da EP terminar no meio de um (falha silenciosa).
  //! lead_in_ms e o trecho antes do primeiro batimento (fase 0).
  static int limit_beats(const TimeBase & tb, int beats,
                         double period_ms, double lead_in_ms = 0.0);

  //! "[t = ... ms | step i/n | tail]"
  static void print_status(const StepContext & s, const std::string & tail);

  //! Abre o arquivo de historico; header vazio = sem cabecalho
  void open_history(const std::string & path, const std::string & header);

  RunPlan           plan_;
  CavityMechanics * mech_ = nullptr;
  std::ofstream     history_;
};

#endif