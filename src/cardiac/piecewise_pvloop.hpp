#ifndef PIECEWISE_PVLOOP_HPP
#define PIECEWISE_PVLOOP_HPP

#include <string>

/*!
 *  ---------------------------------------------------------------------
 *  Controlador de pressao ventricular POR FASES (piecewise PV loop)
 *  ---------------------------------------------------------------------
 *
 *  Implementa o esquema de cinco fases descrito em
 *
 *    Wang, Holmes, Doste, Camps, Margara, Vazquez, Rodriguez (2025/2026),
 *    "Calibration and validation strategy for electromechanical cardiac
 *     digital twins", eLife 14:RP106555, secao "Mechanical model",
 *
 *  que e a formulacao usada no Alya. A diferenca conceitual em relacao ao
 *  caminho -circ (Regazzoni 2022, ver regazzoni2020.hpp) e QUEM e a
 *  incognita:
 *
 *    -circ : o modelo 0D e dono do VOLUME e pergunta ao 3D "que pressao
 *            reproduz este volume?". Isso exige inverter numericamente o
 *            mapa p -> V, ou seja, um Newton 2x2 com varias solucoes 3D
 *            por passo de tempo.
 *
 *    -pw   : o modelo 0D e dono da PRESSAO. A pressao endocardica do passo
 *            n+1 sai de uma equacao explicita que so usa o volume e a taxa
 *            de volume ja medidos no passo n. Uma UNICA solucao 3D por
 *            passo de tempo.
 *
 *  O preco e que as fases isovolumetricas passam a ser impostas por
 *  penalidade (o volume e mantido constante por um controlador PD sobre a
 *  pressao) em vez de por restricao exata, e que a circulacao e aberta:
 *  um Windkessel de dois elementos por ventriculo, sem retorno venoso.
 *
 *  ---------------------------------------------------------------------
 *  UNIDADES INTERNAS: mL, kPa, ms
 *  ---------------------------------------------------------------------
 *  Escolhidas para bater com a Tabela 1 do artigo, onde R_LV aparece em
 *  kPa*ms/mL, C_LV em mL/kPa, as pressoes em kPa e os tempos em ms.
 *
 *      pressao       kPa
 *      volume        mL
 *      tempo         ms
 *      resistencia   kPa*ms/mL
 *      complacencia  mL/kPa
 *
 *  A conversao para a unidade de pressao da mecanica (Pa ou kPa, conforme
 *  -kpa2p) e feita por quem chama, em cardiac_cycle.cpp. Esta classe NAO
 *  conhece a malha nem o solver.
 *
 *  ---------------------------------------------------------------------
 *  AS CINCO FASES (mais a fase 0 de inicializacao)
 *  ---------------------------------------------------------------------
 *
 *  Fase 0 - INIT (inicializacao, so no primeiro batimento)
 *      Rampa linear de 0 ate a pressao de diastase P_DS, ao longo de
 *      t_init. Serve para inflar a malha a partir da configuracao de
 *      referencia sem choque de carga.
 *
 *  Fase 1 - FILLING (enchimento diastolico ativo, "chute atrial")
 *      Rampa linear de P_DS ate EDP ao longo de t_EDP, imitando a
 *      contracao atrial. O fim desta fase e o FIM DA DIASTOLE: o volume
 *      medido ali e V_ED, e e nesse instante que a ativacao eletrica deve
 *      disparar (ver activation_time()).
 *
 *  Fase 2 - IVC (contracao isovolumetrica)
 *      dP = -(1/Cp_iso)(V - V_lock) - (1/Cv_iso) dV/dt
 *      com V_lock = V_ED. O primeiro termo e proporcional ao desvio de
 *      volume, o segundo amortece a taxa; juntos seguram o volume
 *      constante enquanto a tensao ativa sobe.
 *      Saida: quando P alcanca o limiar de ejecao.
 *
 *  Fase 3 - EJECTION (ejecao)
 *      Windkessel de dois elementos:
 *          C_art dP_art/dt + P_art/R_art = -dV/dt,   P = P_art
 *      Saida: quando o fluxo reverte, isto e, dV/dt > 0.
 *
 *      ATENCAO: integrada com ESTABILIZACAO. Tomada ao pe da letra, com
 *      dV/dt medido no passo anterior, esta equacao e instavel: a cavidade
 *      responde a P com V ~ C_cav P, e o modo de periodo 2 do esquema
 *      escalonado tem ganho C_cav/C_art. Para os proprios numeros do
 *      artigo (C_art = 1.39 mL/kPa, elastancia sistolica ~0.45 kPa/mL,
 *      logo C_cav ~ 2.2 mL/kPa) o ganho ja passa de 1.
 *      A cura, na mesma linha do termo de estabilizacao de Regazzoni &
 *      Quarteroni para a tensao ativa, e separar de dV/dt a parte que vem
 *      da propria variacao de pressao e trata-la implicitamente:
 *          (C_art + C_stab) dP/dt + P/R = -[ dV/dt - C_stab dP/dt|_ant ]
 *      O ponto fixo NAO muda (quando dP/dt -> 0 os termos se cancelam e
 *      sobra P/R = -dV/dt); muda so o transiente. Qualquer
 *      C_stab >= C_cav torna o esquema incondicionalmente estavel, e a
 *      complacencia DIASTOLICA medida serve porque superestima a
 *      sistolica -- erra para o lado seguro.
 *
 *  Fase 4 - IVR (relaxamento isovolumetrico)
 *      Mesma equacao da fase 2, com V_lock = V_ES.
 *      Saida: quando P cai abaixo da pressao atrial.
 *
 *  Fase 5 - PASSIVE (enchimento passivo)
 *      dP = -(1/Cp_fill)(V - DSV) - (1/Cv_fill) dV/dt
 *      que devolve o volume ao valor de diastase DSV.
 *      Saida: pelo relogio, no inicio da fase 1 do batimento seguinte.
 *
 *  ---------------------------------------------------------------------
 *  AMBIGUIDADE DE NOTACAO NO ARTIGO, E O QUE FOI FEITO AQUI
 *  ---------------------------------------------------------------------
 *  O artigo escreve as fases 2, 4 e 5 como "dP = ..." (um INCREMENTO de
 *  pressao) e a fase 3 como "C dP/dt + ... = ..." (uma EQUACAO
 *  DIFERENCIAL). Lido ao pe da letra, o ganho das fases isovolumetricas
 *  depende de dt, o que e desconfortavel.
 *
 *  Como nao da para decidir isso pelo texto, as duas leituras estao
 *  implementadas e selecionadas por set_rate_form():
 *
 *      false (default, literal)  P += dP
 *      true  (forma de taxa)     P += dt * dP
 *
 *  No modo de taxa as unidades de Cp e Cv mudam por um fator de tempo
 *  (Cp passa a mL*ms/kPa, Cv a mL/kPa); o log de describe() imprime a
 *  convencao ativa para evitar calibrar contra a errada. A fase 3 e
 *  sempre integrada como equacao diferencial, nas duas leituras.
 *
 *  ---------------------------------------------------------------------
 *  GANHO AUTOMATICO DAS PENALIDADES (Cp/Cv <= 0), E POR QUE E O DEFAULT
 *  ---------------------------------------------------------------------
 *  O artigo NAO tabela Cp e Cv: eles aparecem so como "inverse of penalty
 *  terms". Qualquer valor fixo escolhido aqui seria um chute, e o chute
 *  tem consequencia dura. Linearizando a cavidade em torno do ponto de
 *  travamento, V - V_lock ~ C_cav (P - P_lock), a fase isovolumetrica na
 *  forma literal itera o erro de volume como
 *
 *      e_{n+1} = e_n (1 - C_cav / Cp),
 *
 *  ou seja, ela SO converge se 0 < C_cav/Cp < 2, e e exata (deadbeat) se
 *  Cp = C_cav. C_cav e a complacencia da cavidade, e para os proprios
 *  numeros do artigo ela vale cerca de (105 - 80) mL / (1.3 - 0.6) kPa
 *  ~ 36 mL/kPa. Um Cp fixo da ordem de 20 ja estaria do lado divergente.
 *
 *  Por isso o default e Cp <= 0 = AUTOMATICO: o controlador MEDE a
 *  complacencia diastolica da propria malha (fim da fase 0 contra fim da
 *  fase 1) e usa
 *
 *      dP = -[ (V - V_lock) + damping * dt * dV/dt ] / C_cav,
 *
 *  que e dimensionalmente limpa (mL / (mL/kPa) = kPa), nao depende da
 *  leitura literal/taxa e e sempre estavel: durante a sistole a cavidade
 *  fica MAIS rigida que na diastole, entao C_cav medido na diastole
 *  superestima a complacencia e o ganho erra para o lado lento.
 *
 *  Dar um valor positivo a Cp_iso/Cv_iso/Cp_fill/Cv_fill desliga o modo
 *  automatico daquela penalidade e volta a formula literal do artigo.
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
    double P_DS;        //!< pressao de diastase, alvo da fase 0 [kPa]
    double EDP;         //!< pressao de fim de diastole, alvo da fase 1 [kPa]
    double P_ejection;  //!< limiar de abertura da valva arterial [kPa]
    double P_atrium;    //!< pressao atrial: encerra a fase 4 [kPa]

    //! Volume de diastase alvo da fase 5 [mL]. Negativo = automatico:
    //! usa o volume medido no fim da fase 0, que e por construcao o
    //! estado de diastase desta malha com esta pressao P_DS.
    double DSV;

    double C_art;       //!< complacencia arterial [mL/kPa]
    double R_art;       //!< resistencia arterial [kPa*ms/mL]

    //! Complacencia de ESTABILIZACAO da fase de ejecao [mL/kPa].
    //! <= 0 (default) = usa a complacencia medida C_cav. Ver a nota sobre
    //! a fase 3 no cabecalho da classe. Nao altera o estado estacionario,
    //! so o transiente; existe para tornar o esquema explicito estavel.
    double C_stab;

    //! Penalidades das fases isovolumetricas (2 e 4).
    //! <= 0 (default) = AUTOMATICO a partir da complacencia medida; ver a
    //! nota no cabecalho da classe.
    //! Positivo, forma literal: Cp_iso [mL/kPa], Cv_iso [mL/(kPa*ms)].
    //! Positivo, forma de taxa: Cp_iso [mL*ms/kPa], Cv_iso [mL/kPa].
    double Cp_iso;
    double Cv_iso;

    //! Idem para o enchimento passivo (fase 5).
    double Cp_fill;
    double Cv_fill;

    //! Peso do termo derivativo no modo automatico, adimensional. 0 e o
    //! deadbeat puro; valores maiores amortecem mais e respondem mais
    //! devagar. So tem efeito onde Cp correspondente e <= 0.
    double damping;

    //! Teto de |dP| por passo [kPa]. 0 desativa. E uma rede de seguranca
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
    double P;         //!< pressao endocardica corrente [kPa]
    double P_art;     //!< pressao arterial do Windkessel [kPa]
    double P_art_prev;//!< P_art do passo anterior [kPa]
    double V;         //!< ultimo volume medido no 3D [mL]
    double V_prev;    //!< volume medido no passo anterior [mL]
    double Vdot;      //!< dV/dt medido [mL/ms]
    double V_lock;    //!< volume travado nas fases 2 e 4 [mL]
    double V_DS;      //!< volume de diastase efetivo [mL]
    double V_ED;      //!< volume de fim de diastole [mL]
    double V_ES;      //!< volume de fim de sistole [mL]
    double V_ref;     //!< volume da malha na configuracao de referencia [mL]

    //! Complacencia diastolica MEDIDA da cavidade [mL/kPa]. Estimada no
    //! fim da fase 0 e refinada no fim de cada fase 1. E o ganho do modo
    //! automatico das penalidades; ver o cabecalho da classe.
    double C_cav;
    double t_phase;   //!< instante de entrada na fase corrente [ms]
    int    n_clipped; //!< quantas vezes dP_max atuou (diagnostico)
  };

  PiecewisePVLoop();

  // ---- configuracao -------------------------------------------------

  //! Le o bloco <pvloop_piecewise> da malha (se existir) e em seguida
  //! aplica as sobrescritas de linha de comando. Uma unica funcao de
  //! parsing, de proposito: parametro lido em dois lugares e parametro
  //! que uma hora diverge.
  //!
  //! Todo atributo ausente mantem o default da Tabela 1 do artigo.
  void configure(const std::string & mshfile);

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

  double pressure_LV() const { return st_lv.P; }   //!< [kPa]
  double pressure_RV() const { return st_rv.P; }   //!< [kPa]

  const ChamberState & state_LV() const { return st_lv; }
  const ChamberState & state_RV() const { return st_rv; }

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
                    double V_target, double Cp, double Cv, double dt) const;

  //! Complacencia diastolica a partir de dois pontos (V,P) da propria
  //! malha. Devolve 0 se os pontos nao forem utilizaveis.
  static double compliance_from(double V_hi, double P_hi,
                                double V_lo, double P_lo);
};

#endif
