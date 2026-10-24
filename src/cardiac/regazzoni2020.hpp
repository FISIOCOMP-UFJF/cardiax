#ifndef REGAZZONI2020_HPP
#define REGAZZONI2020_HPP

#include <array>
#include <cmath>
#include <functional>
#include <string>

/*!
 *  Circulacao 0D em malha fechada (Regazzoni et al. 2022): quatro camaras
 *  (LA, LV, RA, RV), quatro valvas e dois circuitos (sistemico e pulmonar),
 *  cada um com compartimentos arterial e venoso RLC.
 *
 *  UNIDADES (as mesmas do resto do cardiax -- nao ha conversao em lugar nenhum):
 *      tempo ms | pressao Pa | volume mL | fluxo mL/ms
 *      resistencia Pa*ms/mL | complacencia mL/Pa | inertancia Pa*ms^2/mL
 *  Os valores publicados (mmHg, s) foram convertidos UMA vez, nos defaults.
 *
 *  Acoplamento 3D-0D: por default cada camara usa elastancia variavel
 *  p = E(t)(V - V0). Com set_BiV_coupling(), as pressoes do VE e do VD vem do
 *  modelo 3D: o 0D entrega os volumes e recebe as pressoes que os reproduzem.
 *
 *  Tempos das camaras: tC/TC/TR sao DERIVADOS (apply_timing) de uma
 *  especificacao invariante com o periodo RR:
 *    - sistole ventricular escala com (RR/RR_ref)^sys_exponent (Bazett, 0.5);
 *    - sistole atrial independe do periodo;
 *    - a contracao atrial e ancorada pelo FIM: termina quando a pressao
 *      ventricular comeca a subir (t_act + emd), o que da sincronia AV correta
 *      em qualquer RR. A variacao do periodo e absorvida pela diastole.
 */
class Regazzoni2020
{
public:
  static const int NSTATE = 12;

  enum StateIdx {
    V_LA = 0, V_LV, V_RA, V_RV,
    p_AR_SYS, p_VEN_SYS, p_AR_PUL, p_VEN_PUL,
    Q_AR_SYS, Q_VEN_SYS, Q_AR_PUL, Q_VEN_PUL
  };

  //! Elastancia variavel de uma camara
  struct Chamber
  {
    double EA;   //!< elastancia ativa  [Pa/mL]
    double EB;   //!< elastancia passiva [Pa/mL]
    double TC;   //!< duracao da contracao [ms]
    double TR;   //!< duracao do relaxamento [ms]
    double tC;   //!< inicio da contracao [ms]
    double V0;   //!< volume em repouso [mL]
  };

  //! Valva tipo diodo: resistencia baixa aberta, alta fechada
  struct Valve { double Rmin, Rmax; };   //!< [Pa*ms/mL]

  //! Circuito arterial + venoso
  struct Compartment
  {
    double R_AR, C_AR, L_AR;
    double R_VEN, C_VEN, L_VEN;
  };

  struct Parameters
  {
    double RR;                    //!< periodo cardiaco [ms]
    Chamber LA, LV, RA, RV;
    Valve MV, AV, TV, PV;
    Compartment SYS, PUL;
  };

  //! Especificacao dos tempos, invariante com o periodo
  struct Timing
  {
    double t_act;            //!< ativacao eletrica ventricular no ciclo [ms]
    double emd;              //!< atraso eletromecanico estimulo -> Ta [ms]
    double TC_atria;         //!< contracao atrial [ms]
    double TR_atria;         //!< relaxamento atrial [ms]
    double TC_vent_ref;      //!< contracao ventricular em RR_ref [ms]
    double TR_vent_ref;      //!< relaxamento ventricular em RR_ref [ms]
    double RR_ref;           //!< periodo a que os dois acima se referem [ms]
    double sys_exponent;     //!< 0.5 = Bazett, 0 = fixa, 1 = proporcional
    double atrial_max_fill;  //!< fracao maxima da diastole que o atrio ocupa
  };

  //! Dados os volumes-alvo e o tempo [ms], devolve as pressoes [Pa] do 3D
  typedef std::function<void(double V_lv, double V_rv, double t,
                             double & p_lv, double & p_rv)> BiVCoupling;

  Regazzoni2020();

  static Parameters default_parameters();
  static Timing     default_timing();
  static std::array<double,NSTATE> default_initial_state();

  const Timing & timing() const { return tim; }
  void set_timing(const Timing & t) { tim = t; apply_timing(); }

  //! Periodo cardiaco [ms]; regera os tempos das camaras.
  void set_period(double rr_ms);
  void set_activation_time(double t_act_ms);   //!< regera os tempos
  void set_emd(double emd_ms);                 //!< regera os tempos

  //! true se a sistole atrial teve de ser encurtada para caber no ciclo
  bool   timing_was_clamped() const { return clamped; }
  double diastolic_window() const;             //!< [ms]; <= 0: ciclo sobrecarregado
  std::string describe_timing() const;

  void set_BiV_coupling(BiVCoupling f) { coupling = f; has_coupling = true; }

  //! Sobrescreve V_LV/V_RV (p. ex. da malha). MUDA o volume total: chame
  //! rebalance_total_volume() depois para preservar a pre-carga.
  void set_initial_volumes(double v_lv, double v_rv)
  { y[V_LV] = v_lv; y[V_RV] = v_rv; }

  //! Devolve o volume total a vtot deslocando por igual a pressao dos dois
  //! reservatorios venosos (os vasos de capacitancia). Retorna o dV [mL].
  double rebalance_total_volume(double vtot);

  double & operator[](int i) { return y[i]; }

  //! Um passo de Euler explicito (condicionalmente estavel: dt ~ 1 ms)
  void step(double t, double dt);

  double period() const { return par.RR; }     //!< [ms]
  double total_volume() const;                 //!< [mL]

  // Grandezas derivadas, validas apos o ultimo rhs()
  double pressure_LA() const { return p_la; }  //!< [Pa]
  double pressure_LV() const { return p_lv; }
  double pressure_RA() const { return p_ra; }
  double pressure_RV() const { return p_rv; }
  double flow_MV() const { return q_mv; }      //!< [mL/ms]
  double flow_AV() const { return q_av; }
  double flow_TV() const { return q_tv; }
  double flow_PV() const { return q_pv; }

  static std::string history_header();
  std::string history_row(double t) const;

private:
  Parameters par;
  Timing     tim;
  std::array<double,NSTATE> y;
  bool clamped;
  BiVCoupling coupling;
  bool has_coupling;
  double p_la, p_lv, p_ra, p_rv;
  double q_mv, q_av, q_tv, q_pv;

  void apply_timing();
  void rhs(double t, const std::array<double,NSTATE> & yy,
           std::array<double,NSTATE> & dy);
  static double wrap(double x, double p);
  double activation(const Chamber & c, double t) const;
  double chamber_pressure(const Chamber & c, double V, double t) const;
  double valve_resistance(const Valve & v, double p1, double p2) const;
  double flux_through_valve(const Valve & v, double p1, double p2) const
  { return (p1 - p2) / valve_resistance(v, p1, p2); }
};

#endif