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
//  UNIDADES, em todo o cardiac_cycle e nos tres modelos, sem conversoes:
//      tempo ms | comprimento um | pressao Pa | volume mL
//  (A unica fronteira com outra unidade e o modelo celular: Cells converte o
//  tempo para a unidade nativa de cada modelo, e CellProfile converte Ta.)

//! O lado 3D, como o modelo de circulacao o enxerga: aplica pressoes nas
//! cavidades e devolve volumes. E a UNICA porta entre os dois mundos.
class CavityMechanics
{
public:
  virtual ~CavityMechanics() {}

  //! UMA solucao mecanica com as pressoes dadas [Pa]
  virtual void apply_pressures(double t_ms, double p_lv, double p_rv) = 0;

  virtual double volume_LV() = 0;            //!< [mL]
  virtual double volume_RV() = 0;            //!< [mL]
  virtual double volume_total() = 0;         //!< cavidade total [mL]

  //! false nos passos em que o log detalhado esta suprimido
  virtual bool verbose() const = 0;
};

//! Entradas de tempo lidas por config() e entregues a configure(). Tudo em ms.
struct TimeBase
{
  double ceiling_ms = -1.0;  //!< -t: TETO da rodada; <= 0 se ausente
  double step_ms    = 1.0;   //!< -dt_mech: passo de acoplamento. O legado o
                             //!< ignora: seu passo e a amostragem da curva.
  int    beats      = 1;     //!< -beats
  double period_ms  = -1.0;  //!< -rr: periodo cardiaco; <= 0 = o de cada modelo
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
  double ep_dt_hint_ms     = 0.0;   //!< passo da EP sugerido pelo modelo (0 = o
                                    //!< do modelo celular); usado se -dt faltar
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