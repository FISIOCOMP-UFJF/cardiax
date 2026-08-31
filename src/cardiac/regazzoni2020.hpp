#ifndef REGAZZONI2020_HPP
#define REGAZZONI2020_HPP

#include <array>
#include <string>
#include <vector>
#include <functional>
#include <cmath>

/*!
 *  Closed-loop lumped-parameter (0D) circulation model.
 *
 *  Regazzoni et al. (2022), "A cardiac electromechanical model coupled with a
 *  lumped-parameter model for closed-loop blood circulation".
 *
 *  Four chambers (LA, LV, RA, RV), four valves (MV, AV, TV, PV) and two
 *  circuits (systemic and pulmonary), each split into arterial and venous
 *  RLC compartments.
 *
 *  ---------------------------------------------------------------------
 *  UNITS: mL, mmHg, s  (consistent with the original model)
 *      resistance  mmHg*s/mL
 *      compliance  mL/mmHg
 *      inertance   mmHg*s^2/mL
 *  The coupling with the 3D mechanics converts pressures at the interface.
 *  ---------------------------------------------------------------------
 *
 *  3D-0D coupling
 *  --------------
 *  By default every chamber uses a time-varying elastance
 *      p(t) = E(t) (V - V0).
 *  When a coupling callback is registered with set_BiV_coupling(), the LV and
 *  RV pressures come instead from the 3D finite element model: the 0D model
 *  hands it the current cavity volumes and receives back the pressures that
 *  reproduce them. The remaining chambers stay 0D.
 *
 *  ---------------------------------------------------------------------
 *  TEMPOS DAS CAMARAS E FREQUENCIA CARDIACA
 *  ---------------------------------------------------------------------
 *  Os instantes tC/TC/TR de cada camara NAO sao mais dados absolutos: eles
 *  sao DERIVADOS de uma especificacao invariante com a frequencia (Timing),
 *  toda vez que set_heart_rate() ou set_timing() e chamado.
 *
 *  A razao e que os tres grupos de tempo escalam de forma diferente:
 *
 *    - sistole ventricular: encurta de forma SUBLINEAR com o periodo,
 *      aproximadamente com (RR/RR_ref)^0.5 (Bazett; a regressao de Weissler
 *      para QS2 da algo muito proximo na faixa 40-120 bpm);
 *
 *    - sistole atrial: e praticamente INDEPENDENTE da frequencia, ~0.17 s;
 *
 *    - inicio da contracao atrial: nao e uma fracao fixa do ciclo. O atrio
 *      tem de TERMINAR de contrair no instante em que a pressao ventricular
 *      comeca a subir. Por isso tC atrial e ancorado pelo FIM,
 *          tC_atrio = t_act + EMD - TC_atrio   (mod RR),
 *      onde t_act e o instante da ativacao eletrica ventricular e EMD e o
 *      atraso eletromecanico (estimulo -> inicio de Ta).
 *
 *  Toda a variacao do periodo e absorvida pela DIASTOLE, que e o que
 *  acontece fisiologicamente. Com isso a sincronia AV fica correta em
 *  qualquer HR, sem reajuste manual.
 *
 *  ATENCAO: nao existe mais "reescalar os tC atuais". Os parametros sao
 *  regerados a partir da especificacao; reescalar um conjunto ja inconsistente
 *  apenas propaga o erro.
 */
class Regazzoni2020
{
public:
  //! Number of ODE state variables
  static const int NSTATE = 12;

  //! State indices
  enum StateIdx {
    V_LA = 0, V_LV, V_RA, V_RV,
    p_AR_SYS, p_VEN_SYS, p_AR_PUL, p_VEN_PUL,
    Q_AR_SYS, Q_VEN_SYS, Q_AR_PUL, Q_VEN_PUL
  };

  //! Time-varying elastance parameters of a 0D chamber
  struct Chamber
  {
    double EA;   //!< active (max) elastance   [mmHg/mL]
    double EB;   //!< passive (min) elastance  [mmHg/mL]
    double TC;   //!< contraction duration     [s]
    double TR;   //!< relaxation duration      [s]
    double tC;   //!< contraction onset        [s]
    double V0;   //!< unloaded volume          [mL]
  };

  //! Diode-like valve: low resistance when open, high when closed
  struct Valve
  {
    double Rmin; //!< [mmHg*s/mL]
    double Rmax; //!< [mmHg*s/mL]
  };

  //! Arterial or venous RLC compartment
  struct Compartment
  {
    double R_AR, C_AR, L_AR;
    double R_VEN, C_VEN, L_VEN;
  };

  struct Parameters
  {
    double HR;                    //!< heart rate [Hz]
    Chamber LA, LV, RA, RV;
    Valve MV, AV, TV, PV;
    Compartment SYS, PUL;
  };

  /*!
   *  Especificacao dos tempos, INVARIANTE com a frequencia cardiaca.
   *  Os tC/TC/TR de Parameters sao gerados a partir daqui por apply_timing().
   */
  struct Timing
  {
    //! Instante da ativacao eletrica ventricular dentro do ciclo [s].
    //! No caminho acoplado e o instante em que o estimulo LAT dispara; o
    //! default 0 corresponde ao laco de solve_circulation(), que estimula
    //! em t = 0, T, 2T, ...
    double t_act;

    //! Atraso eletromecanico: intervalo entre o estimulo e o inicio da subida
    //! de Ta [s]. Propriedade do modelo celular, NAO da frequencia. ~0.04 s
    //! para o ToRORd-Land. E o que define onde a sistole atrial termina.
    double emd;

    //! Duracao da contracao e do relaxamento atriais [s]. Tratadas como
    //! independentes da frequencia (ver comentario do cabecalho da classe).
    double TC_atria;
    double TR_atria;

    //! Duracao da contracao e do relaxamento ventriculares NO PERIODO DE
    //! REFERENCIA [s]. Só entram em jogo sem acoplamento 3D: quando o 3D
    //! fornece p_LV/p_RV a elastancia ventricular e ignorada. Ainda assim
    //! valem a pena porque (a) o modo 0D puro e o diagnostico mais barato
    //! para separar problema de tempo de problema de mecanica e (b) a
    //! validacao da janela diastolica usa esses valores.
    double TC_vent_ref;
    double TR_vent_ref;

    //! Periodo ao qual TC_vent_ref/TR_vent_ref se referem [s].
    double RR_ref;

    //! Expoente do escalonamento sistolico. 0.5 = Bazett (default),
    //! 0.0 = sistole de duracao fixa, 1.0 = escalonamento proporcional.
    double sys_exponent;

    //! Fracao maxima do intervalo diastolico que a sistole atrial pode
    //! ocupar antes de TC_atria ser encurtada. Evita que, em frequencia
    //! alta, a contracao atrial invada a sistole ventricular.
    double atrial_max_fill;
  };

  //! Signature of the 3D coupling: given target cavity volumes and time,
  //! return the LV and RV pressures in mmHg.
  typedef std::function<void(double V_lv, double V_rv, double t,
                             double & p_lv, double & p_rv)> BiVCoupling;

  Regazzoni2020();

  //! Default parameter set from Regazzoni et al.
  //! NOTA: o HR default mudou de 1.25 para 1.0 Hz. Os tC/TC/TR publicados
  //! sao consistentes com um periodo de 1.0 s; o 1.25 anterior colocava a
  //! contracao atrial DENTRO da sistole ventricular (0.90 mod 0.8 = 0.10,
  //! identico ao tC do VE). Os tempos aqui ja saem de apply_timing().
  static Parameters default_parameters();

  //! Especificacao de tempos default (ver Timing).
  static Timing default_timing();

  //! Default (already near periodic) initial state
  static std::array<double,NSTATE> default_initial_state();

  //! Substitui os parametros SEM regerar os tempos. Use quando quiser
  //! controlar tC/TC/TR manualmente (compatibilidade com o comportamento
  //! anterior). Para o caminho generico prefira set_heart_rate().
  void set_parameters(const Parameters & p) { par = p; }
  Parameters & parameters() { return par; }

  const Timing & timing() const { return tim; }

  //! Troca a especificacao de tempos e REGERA os tC/TC/TR.
  void set_timing(const Timing & t) { tim = t; apply_timing(); }

  /*!
   *  Define a frequencia cardiaca [Hz] e REGERA todos os tempos das camaras
   *  a partir da especificacao corrente. Este e o ponto de entrada generico:
   *  qualquer HR produz um ciclo com sincronia AV correta.
   *
   *  Deve ser chamado ANTES de qualquer coisa que dependa de period() --
   *  em particular do dimensionamento do data writer em config(), que usa
   *  circ.period() para contar os frames.
   */
  void set_heart_rate(double hr_hz);

  //! Instante da ativacao eletrica ventricular [s]. Regera os tempos.
  void set_activation_time(double t_act);

  //! Atraso eletromecanico estimulo -> inicio de Ta [s]. Regera os tempos.
  void set_emd(double emd);

  //! Recalcula tC/TC/TR das quatro camaras a partir de tim e de par.HR.
  //! Chamado automaticamente por set_heart_rate/set_timing/set_activation_time.
  void apply_timing();

  //! true se a sistole atrial teve de ser encurtada para caber no ciclo,
  //! ou se a janela diastolica ficou nao positiva. Indica HR alto demais
  //! para a duracao sistolica configurada.
  bool timing_was_clamped() const { return clamped; }

  //! Janela diastolica disponivel no ciclo corrente [s]: o intervalo entre o
  //! fim do relaxamento ventricular e o inicio da contracao atrial, mais a
  //! propria sistole atrial. Negativo significa ciclo sobrecarregado.
  double diastolic_window() const;

  //! Descricao legivel da linha do tempo resultante, para o log.
  std::string describe_timing() const;

  void set_state(const std::array<double,NSTATE> & s) { y = s; }
  const std::array<double,NSTATE> & state() const { return y; }
  double & operator[](int i) { return y[i]; }

  //! Register the 3D mechanics coupling for LV and RV.
  void set_BiV_coupling(BiVCoupling f) { coupling = f; has_coupling = true; }
  void clear_BiV_coupling() { has_coupling = false; }

  //! Override only the initial LV/RV volumes (e.g. from the 3D mesh).
  //! ATENCAO: isto MUDA o volume sanguineo total do modelo. Chame
  //! rebalance_total_volume() em seguida se quiser preservar a pre-carga.
  void set_initial_volumes(double v_lv, double v_rv)
  { y[V_LV] = v_lv; y[V_RV] = v_rv; }

  /*!
   *  Ajusta os reservatorios venosos para que total_volume() volte a valer
   *  vtot. Serve para reparar o volume estressado depois de sobrescrever
   *  V_LV/V_RV com os volumes da malha 3D, que em geral nao coincidem com os
   *  do estado inicial publicado.
   *
   *  O deficit e absorvido por um deslocamento COMUM de pressao nos dois
   *  reservatorios venosos (sistemico e pulmonar), que sao os vasos de
   *  capacitancia. Retorna o volume que foi realocado [mL].
   */
  double rebalance_total_volume(double vtot);

  //! Right-hand side of the ODE system. Fills dy and refreshes the derived
  //! quantities (pressures and flows) accessible through the getters below.
  void rhs(double t, const std::array<double,NSTATE> & yy,
           std::array<double,NSTATE> & dy);

  //! Advance one step. Forward Euler mirrors the reference implementation;
  //! it is conditionally stable, so dt must be small (1e-3 s is typical).
  void step(double t, double dt);

  //! Chamber period
  double period() const { return 1.0 / par.HR; }

  // ---- derived quantities, valid after the last rhs() evaluation ----
  double pressure_LA() const { return p_la; }
  double pressure_LV() const { return p_lv; }
  double pressure_RA() const { return p_ra; }
  double pressure_RV() const { return p_rv; }
  double flow_MV() const { return q_mv; }
  double flow_AV() const { return q_av; }
  double flow_TV() const { return q_tv; }
  double flow_PV() const { return q_pv; }

  //! Total blood volume currently held in the model
  double total_volume() const;

  //! Header and one row for a csv-style history file
  static std::string history_header();
  std::string history_row(double t) const;

private:
  Parameters par;
  Timing     tim;
  std::array<double,NSTATE> y;

  bool clamped;                 //!< a sistole atrial teve de ser encurtada

  BiVCoupling coupling;
  bool has_coupling;

  // derived (static) variables refreshed by rhs()
  double p_la, p_lv, p_ra, p_rv;
  double q_mv, q_av, q_tv, q_pv;

  //! x dobrado em [0, p)
  static double wrap(double x, double p);

  //! Blanco time-varying elastance, normalised activation in [0,1]
  double activation(const Chamber & c, double t) const;

  //! p = E(t) (V - V0)
  double chamber_pressure(const Chamber & c, double V, double t) const;

  //! Smooth diode resistance; ~Rmin when p1 > p2 (open), ~Rmax otherwise
  double valve_resistance(const Valve & v, double p1, double p2) const;

  double flux_through_valve(const Valve & v, double p1, double p2) const
  { return (p1 - p2) / valve_resistance(v, p1, p2); }
};

#endif
