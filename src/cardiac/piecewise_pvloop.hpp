#ifndef PIECEWISE_PVLOOP_HPP
#define PIECEWISE_PVLOOP_HPP

#include <string>

/*!
 *  Controlador de pressao ventricular POR FASES (piecewise PV loop).
 *
 *  Esquema de cinco fases de Wang et al. (2025), eLife 14:RP106555, secao
 *  "Mechanical model" -- a formulacao do Alya. Diferenca em relacao ao
 *  caminho -circ (Regazzoni 2022): QUEM e a incognita.
 *
 *    -circ : o 0D e dono do VOLUME; pergunta ao 3D "que pressao reproduz este
 *            volume?" -> Newton 2x2, varias solucoes 3D por passo.
 *    -pw   : o 0D e dono da PRESSAO. A pressao do passo n+1 sai de uma equacao
 *            explicita que so usa V e dV/dt medidos no passo n -> UMA solucao
 *            3D por passo.
 *
 *  O preco: as fases isovolumetricas sao impostas por penalidade (um PD sobre
 *  a pressao segura o volume) e a circulacao e aberta -- um Windkessel de dois
 *  elementos por ventriculo, sem retorno venoso.
 *
 *  UNIDADES: tempo ms | pressao Pa | volume mL | resistencia Pa*ms/mL |
 *            complacencia mL/Pa. As mesmas do resto do cardiax; esta classe
 *            nao conhece a malha nem o solver. (A Tabela 1 do artigo esta em
 *            kPa; os defaults abaixo ja foram convertidos, x1000.)
 *
 *  AS FASES
 *  0 INIT      rampa de 0 ate P_DS em t_init (infla a malha sem choque).
 *  1 FILLING   rampa de P_DS ate EDP em t_EDP ("chute atrial"). O fim desta
 *              fase e o fim da diastole: V_ED e medido aqui e a ativacao
 *              eletrica dispara neste instante (activation_time()).
 *  2 IVC       contracao isovolumetrica: PD sobre a pressao, V travado em V_ED.
 *              Sai quando P >= limiar de ejecao (P_ejection, ou P_art).
 *  3 EJECTION  Windkessel de 3 elementos: C_art dP_art/dt + P_art/R_art = Q,
 *              P = P_art + R_c Q, Q = -dV/dt. R_c = 0 (default) e o de 2
 *              elementos do artigo. Integrado com C_stab (estabilizacao, ver abaixo).
 *              Sai quando o fluxo reverte (dV/dt > 0), ou P <= P_atrium.
 *  4 IVR       relaxamento isovolumetrico: PD com V travado em V_ES.
 *              Sai quando P <= P_atrium.
 *  5 PASSIVE   enchimento passivo: PD devolve o volume a V_DS, ou, com fill_rate > 0,
 *              segue V = V_entrada + fill_rate (t - t_entrada) (vazao prescrita,
 *              como o PrescribedInflow do pulse). Sai pelo relogio, no inicio da
 *              fase 1 do batimento seguinte.
 *
 *  ESTABILIZACAO DA EJECAO. Com dV/dt do passo anterior, a equacao da fase 3
 *  tomada ao pe da letra e instavel (ganho C_cav/C_art > 1 nos numeros do
 *  artigo). A cura e separar de dV/dt a parte que vem da propria variacao de
 *  pressao e trata-la implicitamente:
 *      (C_art + C_stab) dP/dt + P/R = -[ dV/dt - C_stab dP/dt|ant ]
 *  O ponto fixo nao muda; qualquer C_stab >= C_cav estabiliza.
 *
 *  GANHO AUTOMATICO DAS PENALIDADES (Cp/Cv <= 0, o default). O artigo nao
 *  tabela Cp e Cv. Linearizando a cavidade, a fase isovolumetrica literal
 *  itera e_{n+1} = e_n (1 - C_cav/Cp): so converge se 0 < C_cav/Cp < 2. Por
 *  isso o default MEDE a complacencia diastolica da propria malha e usa
 *      dP = -[ (V - V_lock) + damping*dt*dV/dt ] / C_cav
 *  (deadbeat com amortecimento; sempre estavel, pois a cavidade sistolica e
 *  mais rigida que a diastolica). Um valor positivo em Cp/Cv volta a formula
 *  literal do artigo: dP = -(1/Cp)(V - V_lock) - (1/Cv) dV/dt, com P += dP
 *  (rate_form = false) ou P += dt*dP (rate_form = true; as unidades de Cp e Cv
 *  mudam por um fator de tempo). A fase 3 e sempre uma EDO.
 */
class PiecewisePVLoop
{
public:

  //! Fases do ciclo. Os valores numericos batem com a numeracao do artigo
  //! e sao gravados como coluna do arquivo de historico.
  enum Phase
  {
    INIT     = 0,   //!< rampa de inicializacao ate a diastase
    FILLING  = 1,   //!< enchimento diastolico ativo (contracao atrial)
    IVC      = 2,   //!< contracao isovolumetrica
    EJECTION = 3,   //!< ejecao (Windkessel)
    IVR      = 4,   //!< relaxamento isovolumetrico
    PASSIVE  = 5    //!< enchimento passivo
  };

  //! Nome legivel de uma fase, para o log e para o cabecalho do historico.
  static const char * phase_name(Phase p);

  // -------------------------------------------------------------------
  //  Tempos do ciclo, em ms
  // -------------------------------------------------------------------
  struct Timing
  {
    double RR;         //!< periodo cardiaco [ms] (1000 no artigo)
    double t_init;     //!< duracao da fase 0 [ms]
    double t_EDP;      //!< duracao da fase 1 [ms] (150 no artigo)

    //! Tempo minimo dentro de uma fase antes de aceitar uma transicao
    //! DIRIGIDA POR EVENTO [ms]. Sem isso, o ruido de dV/dt logo apos a
    //! abertura da valva dispara a saida da ejecao no primeiro passo,
    //! porque a taxa herdada da fase isovolumetrica anterior e ~0 e
    //! qualquer jitter positivo ja satisfaz o criterio dV/dt > 0.
    double min_phase;
  };

  // -------------------------------------------------------------------
  //  Parametros de uma cavidade
  // -------------------------------------------------------------------
  struct ChamberParams
  {
    double P_DS;        //!< pressao de diastase, alvo da fase 0 [Pa]
    double EDP;         //!< pressao de fim de diastole, alvo da fase 1 [Pa]
    double P_ejection;  //!< limiar de abertura da valva arterial [Pa]
    double P_atrium;    //!< pressao atrial: encerra a fase 4 [Pa]

    //! Volume de diastase alvo da fase 5 [mL]. Negativo = automatico:
    //! usa o volume medido no fim da fase 0, que e por construcao o
    //! estado de diastase desta malha com esta pressao P_DS.
    double DSV;

    double C_art;       //!< complacencia arterial [mL/Pa]
    double R_art;       //!< resistencia arterial [Pa*ms/mL]

    //! Impedancia caracteristica [Pa*ms/mL]: resistencia em serie entre a
    //! cavidade e a complacencia arterial, P = P_art + R_c Q. 0 = Windkessel
    //! de 2 elementos (o do artigo).
    double R_c;

    //! Complacencia de ESTABILIZACAO da fase de ejecao [mL/Pa].
    //! <= 0 (default) = usa a complacencia medida C_cav. Ver a nota sobre
    //! a fase 3 no cabecalho da classe. Nao altera o estado estacionario,
    //! so o transiente; existe para tornar o esquema explicito estavel.
    double C_stab;

    //! Penalidades das fases isovolumetricas (2 e 4).
    //! <= 0 (default) = AUTOMATICO a partir da complacencia medida; ver a
    //! nota no cabecalho da classe.
    //! Positivo, forma literal: Cp_iso [mL/Pa], Cv_iso [mL/(Pa*ms)].
    //! Positivo, forma de taxa: Cp_iso [mL*ms/Pa], Cv_iso [mL/Pa].
    double Cp_iso;
    double Cv_iso;

    //! Idem para o enchimento passivo (fase 5).
    double Cp_fill;
    double Cv_fill;

    //! Peso do termo derivativo no modo automatico, adimensional. 0 e o
    //! deadbeat puro; valores maiores amortecem mais e respondem mais
    //! devagar. So tem efeito onde Cp correspondente e <= 0.
    double damping;

    //! Antecipacao no modo automatico das penalidades (fases 2 e 4): alem
    //! do desvio de volume, compensa a parte de dV/dt que NAO vem da variacao
    //! de pressao do passo anterior (encurtamento ou relaxamento pela Ta),
    //!     dP = -[ (V - V_alvo) + dt (dV/dt - C_cav dP_ant/dt) ] / C_cav.
    //! Substitui o termo de damping. Estavel enquanto C_real/C_cav < ~4/3.
    bool feedforward;

    //! Vazao de enchimento da fase 5 [mL/ms]. <= 0 = PD de volta a V_DS.
    double fill_rate;

    //! Teto de |dP| por passo [Pa]. 0 desativa. E uma rede de seguranca
    //! numerica, nao um parametro fisiologico: se ela dispara com
    //! frequencia, o ganho da penalidade esta alto demais para o dt.
    double dP_max;
  };

  // -------------------------------------------------------------------
  //  Estado de uma cavidade
  // -------------------------------------------------------------------
  struct ChamberState
  {
    Phase  phase;
    double P;         //!< pressao endocardica corrente [Pa]
    double P_art;     //!< pressao arterial do Windkessel [Pa]
    double P_prev;    //!< pressao da cavidade no passo anterior [Pa]
    double V;         //!< ultimo volume medido no 3D [mL]
    double V_prev;    //!< volume medido no passo anterior [mL]
    double Vdot;      //!< dV/dt medido [mL/ms]
    double V_lock;    //!< volume travado nas fases 2 e 4; volume de entrada na 5 [mL]
    double V_DS;      //!< volume de diastase efetivo [mL]
    double V_ED;      //!< volume de fim de diastole [mL]
    double V_ES;      //!< volume de fim de sistole [mL]
    double V_ref;     //!< volume da malha na configuracao de referencia [mL]
    bool   has_ejected; //!< ja houve ao menos uma ejecao (-pwwk 2)

    //! Complacencia diastolica MEDIDA da cavidade [mL/Pa]. Estimada no
    //! fim da fase 0 e refinada no fim de cada fase 1. E o ganho do modo
    //! automatico das penalidades; ver o cabecalho da classe.
    double C_cav;
    double t_phase;   //!< instante de entrada na fase corrente [ms]
    int    n_clipped; //!< quantas vezes dP_max atuou (diagnostico)
  };

  PiecewisePVLoop();

  // ---- configuracao -------------------------------------------------

  //! Le o bloco <pvloop_piecewise> da malha (se existir) e em seguida as
  //! flags -pw*. Atributo ausente = default da Tabela 1 do artigo.
  //! rr_ms > 0 (a flag global -rr) sobrescreve o periodo do XML.
  void configure(const std::string & mshfile, double rr_ms = -1.0);

  Timing       & timing() { return tim; }
  ChamberParams & lv()    { return par_lv; }
  ChamberParams & rv()    { return par_rv; }

  const Timing       & timing() const { return tim; }
  const ChamberParams & lv()    const { return par_lv; }
  const ChamberParams & rv()    const { return par_rv; }

  //! Ver a nota sobre a ambiguidade de notacao no cabecalho da classe.
  void set_rate_form(bool r) { rate_form = r; }
  bool rate_form_active() const { return rate_form; }

  //! Windkessel continuo (true) contra limiar fixo (false, default).
  //!
  //! false: a fase 3 comeca quando P >= P_ejection, e P_art e reiniciado
  //!        em P_ejection na entrada da fase. E a leitura literal do
  //!        artigo, em que "ejection pressure" e um parametro calibrado.
  //!
  //! true : P_art e um estado persistente, integrado tambem fora da
  //!        ejecao com C dP_art/dt = -P_art/R_art (runoff diastolico), e
  //!        a fase 3 comeca quando P >= P_art. E o comportamento
  //!        fisiologico da aorta, e faz a pressao diastolica arterial
  //!        emergir do modelo em vez de ser prescrita.
  //!
  //! hold_until_ejection (-pwwk 2): como true, mas P_art so comeca a drenar
  //!        depois da PRIMEIRA ejecao da cavidade (P_art = P_ejection ate la),
  //!        como o pulse.cycle (has_ejected).
  void set_continuous_windkessel(bool c) { continuous_wk = c; }

  // ---- ciclo de vida ------------------------------------------------

  //! Estado inicial. V_lv0/V_rv0 sao os volumes da malha na configuracao
  //! de referencia, medidos ANTES de qualquer solucao mecanica.
  void init(double V_lv0, double V_rv0);

  //! Calcula as pressoes do passo que TERMINA em t_end, a partir do
  //! volume e da taxa medidos no passo anterior. Chamar ANTES da
  //! solucao 3D.
  void advance_pressure(double t_end, double dt);

  //! Registra os volumes medidos no 3D no instante t e avalia as
  //! transicoes de fase. Chamar DEPOIS da solucao 3D.
  void update_volumes(double t, double dt, double V_lv, double V_rv);

  // ---- consultas ----------------------------------------------------

  double pressure_LV() const { return st_lv.P; }   //!< [Pa]
  double pressure_RV() const { return st_rv.P; }   //!< [Pa]

  const ChamberState & state_LV() const { return st_lv; }
  const ChamberState & state_RV() const { return st_rv; }

  //! Repoe o estado das cavidades (restauracao de checkpoint)
  void set_states(const ChamberState & lv, const ChamberState & rv)
  { st_lv = lv; st_rv = rv; }

  //! Instante da ativacao eletrica do batimento b (0-based), em ms:
  //! o fim da fase 1, que e onde o artigo posiciona o inicio da fase 2.
  double activation_time(int beat) const
  { return tim.t_init + beat * tim.RR + tim.t_EDP; }

  //! Deslocamento a aplicar ao relogio do estimulo LAT, em ms. E o
  //! instante da PRIMEIRA ativacao; as demais saem dele por periodicidade.
  double stimulus_offset() const { return activation_time(0); }

  //! Duracao total da rodada com num_beats batimentos, em ms.
  double total_time(int num_beats) const
  { return tim.t_init + num_beats * tim.RR; }

  //! Batimento a que pertence o instante t (0-based). Negativo durante a
  //! fase de inicializacao.
  int beat_of(double t) const;

  //! Volume sistolico e fracao de ejecao do ultimo ciclo completo do VE.
  double stroke_volume_LV() const { return st_lv.V_ED - st_lv.V_ES; }
  double ejection_fraction_LV() const;
  double stroke_volume_RV() const { return st_rv.V_ED - st_rv.V_ES; }
  double ejection_fraction_RV() const;

  //! Resumo legivel da configuracao, para o log de config().
  std::string describe() const;

  //! Cabecalho e uma linha do arquivo de historico.
  static std::string history_header();
  std::string history_row(double t) const;

private:
  Timing        tim;
  ChamberParams par_lv, par_rv;
  ChamberState  st_lv,  st_rv;

  bool rate_form;
  bool continuous_wk;
  bool wk_hold_until_ejection;

  //! Instante em que comeca a fase 1 do batimento b [ms].
  double filling_start(int beat) const
  { return tim.t_init + beat * tim.RR; }

  //! Aplica o teto dP_max e contabiliza.
  double clip_dP(double dP, ChamberState & st,
                 const ChamberParams & p) const;

  //! Pressao da cavidade no passo corrente. Nao toca no volume.
  void advance_one(ChamberState & st, const ChamberParams & p,
                   double t_end, double dt) const;

  //! Transicoes de fase de uma cavidade.
  void transition_one(ChamberState & st, const ChamberParams & p,
                      double t, const char * tag);

  //! Runoff diastolico do Windkessel, so no modo continuo.
  void arterial_runoff(ChamberState & st, const ChamberParams & p,
                       double dt) const;

  //! Incremento de pressao de uma fase controlada por penalidade (2, 4 ou
  //! 5). Escolhe entre o modo automatico (Cp <= 0) e a formula literal do
  //! artigo, e ja aplica o teto dP_max.
  double penalty_dP(ChamberState & st, const ChamberParams & p,
                    double V_target, double Cp, double Cv, double dt,
                    bool allow_ff = true) const;

  //! Complacencia diastolica a partir de dois pontos (V,P) da propria
  //! malha. Devolve 0 se os pontos nao forem utilizaveis.
  static double compliance_from(double V_hi, double P_hi,
                                double V_lo, double P_lo);
};

#endif