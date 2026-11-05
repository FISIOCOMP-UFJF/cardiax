#ifndef PIECEWISE_PVLOOP_HPP
#define PIECEWISE_PVLOOP_HPP

#include <string>

/*!
 *  Ciclo cardiaco de cinco fases, estilo Alya: PORT DIRETO do pulse.cycle
 *  (fenicsx-pulse 0.12, pulse/cycle.py), sem ajustes proprios do cardiax.
 *
 *  Origem do esquema: Wang et al. (2025), eLife 14:RP106555, secao
 *  "Mechanical model". A implementacao seguida aqui, linha a linha, e a do
 *  pulse (que por sua vez e um port do cycle_controller do physcardems).
 *
 *  Cada cavidade tem UMA incognita de pressao. Em cada passo, a fase corrente
 *  impoe UMA restricao sobre ela, calculada com o estado do INICIO do passo:
 *
 *    0 PRELOAD    pressao: rampa linear (0 -> preload_pressure em t_zero no
 *                 1o batimento; depois preload_pressure -> p_end_diastole
 *                 entre t_zero e t_end_diastole de cada batimento)
 *    1 IVC        volume:  V = V_ED
 *    2 EJECTION   pressao afim em V: Windkessel de 3 elementos, Euler
 *                 implicito,
 *                     Q   = -(V - V_n)/dt
 *                     P_c = (P_c,n + dt Q/C) / (1 + dt/(R_p C))
 *                     P_v = P_c + R_c Q          ->  P_v = A + B V
 *    3 IVR        volume:  V = V_ES
 *    4 FILLING    volume:  V = V_n + rate*dt (vazao prescrita)
 *
 *  Transicoes, avaliadas com o estado do FIM do passo (ja convergido):
 *    PRELOAD -> IVC      t >= n_beats*periodo + t_end_diastole
 *    IVC -> EJECTION     P > P_c (abre a valva arterial)
 *    EJECTION -> IVR     passados mais de min_ejection_duration desde a
 *                        abertura: dV > -dvol_eps (fluxo parou) ou P < P_c
 *    IVR -> FILLING      P < p_fill
 *    FILLING -> PRELOAD  t >= (n_beats+1)*periodo + t_zero
 *  Com a valva fechada, P_c fica parado ate a primeira ejecao da cavidade e
 *  depois drena pela resistencia periferica (Euler implicito).
 *
 *  Esta classe e so a POLITICA: nao conhece a malha nem o solver. Quem impoe a
 *  restricao (pulse: monoliticamente no Newton; cardiax: Newton externo nas
 *  pressoes, ver PiecewiseCirculation) e quem chama constraint() e commit().
 *
 *  UNIDADES: tempo ms | pressao Pa | volume mL | resistencia Pa*ms/mL |
 *            complacencia mL/Pa | vazao mL/ms. (O pulse usa SI: s, Pa, m^3.)
 */
class PiecewisePVLoop
{
public:

  //! Fases, com a numeracao do pulse.cycle.Phase
  enum Phase
  {
    PRELOAD  = 0,
    IVC      = 1,   //!< ISOVOLUMIC_CONTRACTION
    EJECTION = 2,
    IVR      = 3,   //!< ISOVOLUMIC_RELAXATION
    FILLING  = 4
  };

  static const char * phase_name(Phase p);

  //! Windkessel de 3 elementos (pulse.cycle.Windkessel). Sem estado: P_c
  //! vive no ChamberState.
  struct Windkessel
  {
    double p_init;   //!< P_c inicial [Pa]
    double C;        //!< complacencia [mL/Pa]
    double R;        //!< resistencia periferica [Pa*ms/mL]
    double Rc;       //!< impedancia caracteristica [Pa*ms/mL] (0 = 2 elementos)

    //! Coeficientes de P_v = A + B V impostos durante a ejecao
    void ejection_law(double P_c, double V_n, double dt, double & A, double & B) const;

    //! Avanca P_c num passo convergido; devolve (P_c novo, Q)
    void advance(double P_c, double V, double V_n, double dt, bool ejecting,
                 bool has_ejected, double & P_c_new, double & Q) const;
  };

  //! pulse.cycle.CycleParams de uma cavidade
  struct ChamberParams
  {
    double t_zero;                 //!< [ms]
    double preload_pressure;       //!< [Pa]
    double t_end_diastole;         //!< [ms], medido do inicio do batimento
    double p_end_diastole;         //!< [Pa]
    double p_fill;                 //!< [Pa]
    double period;                 //!< [ms]
    Windkessel windkessel;
    double fill_rate;              //!< PrescribedInflow.rate [mL/ms]
    double min_ejection_duration;  //!< [ms] (pulse: 0.01 s)
    double dvol_eps;               //!< [mL]  (pulse: 1e-10 m^3 = 1e-4 mL)
  };

  //! pulse.cycle.CavityCycle, mais a fase sob a qual o ultimo passo foi resolvido
  struct ChamberState
  {
    Phase  phase;
    int    n_beats;
    double last_phase_change;   //!< [ms]
    double volume_n;            //!< [mL]
    double pressure_n;          //!< [Pa]
    double dvol_n;              //!< [mL]
    double wdk_pressure_n;      //!< P_c [Pa]
    double outflow_n;           //!< Q [mL/ms]
    bool   has_ejected;
    double end_dia_vol;         //!< [mL]
    double end_sys_vol;         //!< [mL]
    Phase  solved_under;        //!< fase do ultimo passo resolvido (historico)
  };

  //! Restricao de um passo sobre a pressao de uma cavidade (CavityControl)
  struct Constraint
  {
    enum Kind { PRESSURE, VOLUME, AFFINE } kind;
    double P;      //!< PRESSURE: p = P
    double V;      //!< VOLUME:   V(p) = V
    double A, B;   //!< AFFINE:   p = A + B V(p)   [Pa, Pa/mL]
  };

  PiecewisePVLoop();

  // ---- configuracao -------------------------------------------------

  //! Bloco <pvloop_piecewise> da malha (se existir), depois as flags -pw*.
  //! rr_ms > 0 (flag global -rr) sobrescreve o periodo.
  void configure(const std::string & mshfile, double rr_ms = -1.0);

  ChamberParams & lv() { return par_lv; }
  ChamberParams & rv() { return par_rv; }
  const ChamberParams & lv() const { return par_lv; }
  const ChamberParams & rv() const { return par_rv; }

  // ---- ciclo de vida (CycleController) ------------------------------

  //! initialize(t0): volume e pressao correntes, todas as cavidades em PRELOAD
  void initialize(double t0, double V_lv, double P_lv, double V_rv, double P_rv);

  //! Restricoes do passo que TERMINA em t, a partir do estado do inicio dele
  Constraint constraint_LV(double t, double dt) const { return constraint(st_lv, par_lv, t, dt); }
  Constraint constraint_RV(double t, double dt) const { return constraint(st_rv, par_rv, t, dt); }

  //! Passo convergido em t: avanca o Windkessel e decide a proxima fase.
  //! Devolve true se alguma cavidade mudou de fase.
  bool commit(double t, double dt, double V_lv, double P_lv, double V_rv, double P_rv);

  // ---- consultas ----------------------------------------------------

  const ChamberState & state_LV() const { return st_lv; }
  const ChamberState & state_RV() const { return st_rv; }
  void set_states(const ChamberState & lv, const ChamberState & rv) { st_lv = lv; st_rv = rv; }

  double period() const { return par_lv.period; }

  //! Inicio da ativacao do 1o batimento [ms]: o pulse poe t_sys da Bestel em
  //! t_end_diastole.
  double activation_time() const { return par_lv.t_end_diastole; }

  //! Duracao da rodada (o pulse roda NUM_BEATS * PERIOD) [ms]
  double total_time(int num_beats) const { return num_beats * period(); }

  //! Batimento do passo que termina em t (0-based)
  int beat_of(double t) const;

  //! VS/FE com os V_ED/V_ES mais recentes
  double stroke_volume_LV() const { return st_lv.end_dia_vol - st_lv.end_sys_vol; }
  double stroke_volume_RV() const { return st_rv.end_dia_vol - st_rv.end_sys_vol; }
  double ejection_fraction_LV() const;
  double ejection_fraction_RV() const;

  std::string describe() const;

  //! Historico: uma linha por passo, como o history do complete_cycle.py
  //! (a fase gravada e a fase sob a qual o passo foi RESOLVIDO).
  static std::string history_header();
  std::string history_row(double t) const;

private:
  ChamberParams par_lv, par_rv;
  ChamberState  st_lv,  st_rv;

  static Constraint constraint(const ChamberState & st, const ChamberParams & p,
                               double t, double dt);
  static double preload_pressure(const ChamberState & st, const ChamberParams & p, double t);
  static bool advance_phase(ChamberState & st, const ChamberParams & p, double t,
                            const char * tag);
  static void commit_one(ChamberState & st, const ChamberParams & p, double t, double dt,
                         double V, double P, const char * tag, bool & changed);
};

#endif
