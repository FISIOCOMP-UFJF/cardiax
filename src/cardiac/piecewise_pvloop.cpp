#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "piecewise_pvloop.hpp"
#include "util/command_line_args.h"
#include "util/pugixml.hpp"

// ======================================================================
//  Utilitarios locais
// ======================================================================

namespace
{
  //! Le um atributo XML se ele existir, senao devolve o default. Escrito a
  //! mao em vez de usar as_double(def) porque queremos distinguir "ausente"
  //! de "presente e igual a zero" no relatorio.
  double xml_attr(const pugi::xml_node & node, const char * name, double def)
  {
    if (!node) return def;
    pugi::xml_attribute a = node.attribute(name);
    return a ? a.as_double() : def;
  }

  int xml_attr_i(const pugi::xml_node & node, const char * name, int def)
  {
    if (!node) return def;
    pugi::xml_attribute a = node.attribute(name);
    return a ? a.as_int() : def;
  }

  double clamp01(double x)
  {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
  }
}

// ======================================================================
//  Construcao e defaults
// ======================================================================

PiecewisePVLoop::PiecewisePVLoop()
  : rate_form(false), continuous_wk(false), wk_hold_until_ejection(false)
{
  // ---- tempos, Tabela 1 do artigo ------------------------------------
  tim.RR        = 1000.0;   // "Cycle length 1000 ms"
  tim.t_EDP     =  150.0;   // "Atrial filling duration LV ~150 ms"
  tim.t_init    =  100.0;   // nao tabelado: rampa de inflacao inicial
  tim.min_phase =   10.0;   // guarda numerica, nao fisiologica

  // ---- ventriculo esquerdo -------------------------------------------
  // Valores da coluna "Initial value" da Tabela 1 quando existem. Os
  // valores CALIBRADOS da Figura 1F estao no README e no XML de exemplo;
  // sao eles que produzem o laco fisiologico do artigo, mas os iniciais
  // sao o ponto de partida documentado.
  par_lv.P_DS       =   600.0;   // nao tabelado (diastase ~4.5 mmHg)
  par_lv.EDP        =  1300.0;   // "EDP_LV 1.3 kPa"
  par_lv.P_ejection = 10000.0;   // "P ejection LV 10 kPa" (calibrado: 8450)
  par_lv.P_atrium   =  1300.0;   // "P diastolic cut-off LV 1.3 kPa"
  par_lv.DSV        = -1.0;      // automatico: volume no fim da fase 0
  par_lv.C_art      =  0.0015;   // "C_LV 1.5 mL/kPa" (calibrado: 0.001392)
  par_lv.R_art      = 10000.0;   // "R_LV 10 kPa ms/mL" (calibrado: 21590)
  par_lv.R_c        = 0.0;       // nao existe no artigo (Windkessel de 2 elementos)
  par_lv.C_stab     = -1.0;   // automatico: usa a C_cav medida

  // Penalidades das fases isovolumetricas e do enchimento passivo. NAO
  // sao tabeladas no artigo -- la elas aparecem apenas como "inverse of
  // penalty terms". Negativo = automatico a partir da complacencia medida
  // da propria cavidade; ver a nota no cabecalho de piecewise_pvloop.hpp.
  // Um valor FIXO aqui seria um chute com consequencia dura: a fase
  // isovolumetrica so converge se C_cav/Cp < 2.
  par_lv.Cp_iso     = -1.0;
  par_lv.Cv_iso     = -1.0;
  par_lv.Cp_fill    = -1.0;
  par_lv.Cv_fill    = -1.0;
  par_lv.damping    =  0.3;
  par_lv.feedforward = false;
  par_lv.fill_rate  = -1.0;      // nao tabelado: o artigo devolve a V_DS
  par_lv.dP_max     = 5000.0;

  // ---- ventriculo direito ---------------------------------------------
  par_rv.P_DS       =   400.0;
  par_rv.EDP        =  1000.0;   // "EDP_RV 1.0 kPa"
  par_rv.P_ejection =  2300.0;   // "P ejection RV 2.3 kPa"
  par_rv.P_atrium   =  1500.0;   // "P diastolic cut-off RV 1.5 kPa"
  par_rv.DSV        = -1.0;
  par_rv.C_art      =  0.0045;   // "C_RV 4.5 mL/kPa" (3x o do VE)

  // ATENCAO: a Tabela 1 imprime R_RV = 33.3 kPa ms/mL e descreve o valor como
  // "a third of left ventricular value", o que com R_LV = 10 daria 3.33. O
  // default e 3330 Pa*ms/mL: uma pos-carga pulmonar 3x maior que a sistemica
  // inverteria a relacao de pressoes (o artigo reporta pico de 2.3 kPa no VD
  // contra 14 kPa no VE). Para o numero impresso, R_art="33300" no XML.
  par_rv.R_art      =  3330.0;
  par_rv.R_c        = 0.0;
  par_rv.C_stab     = -1.0;

  par_rv.Cp_iso     = -1.0;
  par_rv.Cv_iso     = -1.0;
  par_rv.Cp_fill    = -1.0;
  par_rv.Cv_fill    = -1.0;
  par_rv.damping    =  0.3;
  par_rv.feedforward = false;
  par_rv.fill_rate  = -1.0;
  par_rv.dP_max     = 5000.0;
}

const char * PiecewisePVLoop::phase_name(Phase p)
{
  switch (p)
  {
    case INIT:     return "INIT";
    case FILLING:  return "FILLING";
    case IVC:      return "IVC";
    case EJECTION: return "EJECTION";
    case IVR:      return "IVR";
    case PASSIVE:  return "PASSIVE";
  }
  return "?";
}

// ======================================================================
//  Configuracao: XML primeiro, linha de comando depois
// ======================================================================

void PiecewisePVLoop::configure(const std::string & mshfile, double rr_ms)
{
  // ---- 1. bloco <pvloop_piecewise> da malha, se existir ---------------
  //
  // Deliberadamente um bloco NOVO, e nao uma extensao do <pvloop> legado:
  // o <pvloop> antigo carrega a curva de pressao tabelada e os parametros
  // do modelo de Windkessel fechado do caminho classico. Misturar os dois
  // conjuntos no mesmo no faria com que uma malha valida para um caminho
  // parecesse valida para o outro.
  pugi::xml_document doc;
  pugi::xml_parse_result res = doc.load_file(mshfile.c_str());

  pugi::xml_node root;
  if (res)
    root = doc.child("mesh").child("pvloop_piecewise");

  if (root)
    std::cout << " <pvloop_piecewise> encontrado na malha" << std::endl;
  else
    std::cout << " <pvloop_piecewise> ausente na malha; usando os defaults"
              << " da Tabela 1 (Wang et al. 2025)" << std::endl;

  tim.RR        = xml_attr(root, "RR",        tim.RR);
  tim.t_init    = xml_attr(root, "t_init",    tim.t_init);
  tim.t_EDP     = xml_attr(root, "t_EDP",     tim.t_EDP);
  tim.min_phase = xml_attr(root, "min_phase", tim.min_phase);

  rate_form     = (xml_attr_i(root, "rate_form",     rate_form     ? 1 : 0) != 0);
  const int wk_xml = xml_attr_i(root, "continuous_wk", continuous_wk ? 1 : 0);
  continuous_wk          = (wk_xml != 0);
  wk_hold_until_ejection = (wk_xml == 2);

  pugi::xml_node nlv = root ? root.child("lv") : pugi::xml_node();
  pugi::xml_node nrv = root ? root.child("rv") : pugi::xml_node();

  ChamberParams * pp[2] = { &par_lv, &par_rv };
  pugi::xml_node  nn[2] = { nlv, nrv };

  for (int c = 0; c < 2; c++)
  {
    ChamberParams & p = *pp[c];
    const pugi::xml_node & n = nn[c];

    p.P_DS       = xml_attr(n, "P_DS",       p.P_DS);
    p.EDP        = xml_attr(n, "EDP",        p.EDP);
    p.P_ejection = xml_attr(n, "P_ejection", p.P_ejection);
    p.P_atrium   = xml_attr(n, "P_atrium",   p.P_atrium);
    p.DSV        = xml_attr(n, "DSV",        p.DSV);
    p.C_art      = xml_attr(n, "C_art",      p.C_art);
    p.R_art      = xml_attr(n, "R_art",      p.R_art);
    p.R_c        = xml_attr(n, "R_c",        p.R_c);
    p.C_stab     = xml_attr(n, "C_stab",     p.C_stab);
    p.Cp_iso     = xml_attr(n, "Cp_iso",     p.Cp_iso);
    p.Cv_iso     = xml_attr(n, "Cv_iso",     p.Cv_iso);
    p.Cp_fill    = xml_attr(n, "Cp_fill",    p.Cp_fill);
    p.Cv_fill    = xml_attr(n, "Cv_fill",    p.Cv_fill);
    p.damping    = xml_attr(n, "damping",    p.damping);
    p.fill_rate  = xml_attr(n, "fill_rate",  p.fill_rate);
    p.feedforward = (xml_attr_i(n, "feedforward", p.feedforward ? 1 : 0) != 0);
    p.dP_max     = xml_attr(n, "dP_max",     p.dP_max);
  }

  // ---- 2. sobrescritas de linha de comando ----------------------------
  //
  // Periodo: -rr (flag global) > XML > default
  if (rr_ms > 0.0) tim.RR = rr_ms;

  tim.t_init    = CommandLineArgs::read("-pwinit", tim.t_init,
      "duracao da fase 0 (inflacao ate a diastase), em ms",
      "pvloop por fases");
  tim.t_EDP     = CommandLineArgs::read("-pwedp", tim.t_EDP,
      "duracao da fase 1 (chute atrial), em ms", "pvloop por fases");
  tim.min_phase = CommandLineArgs::read("-pwminphase", tim.min_phase,
      "tempo minimo em uma fase antes de aceitar transicao por evento, ms",
      "pvloop por fases");

  rate_form = (CommandLineArgs::read("-pwrate", rate_form ? 1 : 0,
      "1 = fases 2/4/5 integradas como dP/dt; 0 = dP literal do artigo",
      "pvloop por fases") != 0);
  const int wk = CommandLineArgs::read("-pwwk", wk_hold_until_ejection ? 2 : (continuous_wk ? 1 : 0),
      "0 = limiar fixo P_ejection; 1 = Windkessel continuo; 2 = continuo, drena"
      " so apos a 1a ejecao (pulse)", "pvloop por fases");
  continuous_wk          = (wk != 0);
  wk_hold_until_ejection = (wk == 2);

  par_lv.P_DS       = CommandLineArgs::read("-pwpds",     par_lv.P_DS,
      "pressao de diastase do VE, Pa", "pvloop por fases");
  par_lv.EDP        = CommandLineArgs::read("-pwedplv",   par_lv.EDP,
      "pressao de fim de diastole do VE, Pa", "pvloop por fases");
  par_lv.P_ejection = CommandLineArgs::read("-pwpej",     par_lv.P_ejection,
      "limiar de ejecao do VE, Pa", "pvloop por fases");
  par_lv.P_atrium   = CommandLineArgs::read("-pwpla",     par_lv.P_atrium,
      "pressao atrial esquerda, Pa", "pvloop por fases");
  par_lv.DSV        = CommandLineArgs::read("-pwdsv",     par_lv.DSV,
      "volume de diastase do VE, mL (<0 = medido no fim da fase 0)",
      "pvloop por fases");
  par_lv.C_art      = CommandLineArgs::read("-pwcart",    par_lv.C_art,
      "complacencia aortica, mL/Pa", "pvloop por fases");
  par_lv.R_art      = CommandLineArgs::read("-pwrart",    par_lv.R_art,
      "resistencia aortica, Pa*ms/mL", "pvloop por fases");
  par_lv.R_c        = CommandLineArgs::read("-pwrc",      par_lv.R_c,
      "impedancia caracteristica aortica, Pa*ms/mL (0 = Windkessel de 2 elementos)",
      "pvloop por fases");
  par_lv.Cp_iso     = CommandLineArgs::read("-pwcpiso",   par_lv.Cp_iso,
      "penalidade de volume das fases isovolumetricas do VE",
      "pvloop por fases");
  par_lv.Cv_iso     = CommandLineArgs::read("-pwcviso",   par_lv.Cv_iso,
      "penalidade de taxa das fases isovolumetricas do VE",
      "pvloop por fases");
  par_lv.Cp_fill    = CommandLineArgs::read("-pwcpfill",  par_lv.Cp_fill,
      "penalidade de volume do enchimento passivo do VE",
      "pvloop por fases");
  par_lv.Cv_fill    = CommandLineArgs::read("-pwcvfill",  par_lv.Cv_fill,
      "penalidade de taxa do enchimento passivo do VE", "pvloop por fases");

  par_rv.P_DS       = CommandLineArgs::read("-pwpds_rv",  par_rv.P_DS,
      "pressao de diastase do VD, Pa", "pvloop por fases");
  par_rv.EDP        = CommandLineArgs::read("-pwedprv",   par_rv.EDP,
      "pressao de fim de diastole do VD, Pa", "pvloop por fases");
  par_rv.P_ejection = CommandLineArgs::read("-pwpej_rv",  par_rv.P_ejection,
      "limiar de ejecao do VD, Pa", "pvloop por fases");
  par_rv.P_atrium   = CommandLineArgs::read("-pwpra",     par_rv.P_atrium,
      "pressao atrial direita, Pa", "pvloop por fases");
  par_rv.DSV        = CommandLineArgs::read("-pwdsv_rv",  par_rv.DSV,
      "volume de diastase do VD, mL (<0 = medido no fim da fase 0)",
      "pvloop por fases");
  par_rv.C_art      = CommandLineArgs::read("-pwcart_rv", par_rv.C_art,
      "complacencia pulmonar, mL/Pa", "pvloop por fases");
  par_rv.R_art      = CommandLineArgs::read("-pwrart_rv", par_rv.R_art,
      "resistencia pulmonar, Pa*ms/mL", "pvloop por fases");
  par_rv.R_c        = CommandLineArgs::read("-pwrc_rv",   par_rv.R_c,
      "impedancia caracteristica pulmonar, Pa*ms/mL (0 = Windkessel de 2 elementos)",
      "pvloop por fases");
  par_rv.Cp_iso     = CommandLineArgs::read("-pwcpiso_rv", par_rv.Cp_iso,
      "penalidade de volume das fases isovolumetricas do VD",
      "pvloop por fases");
  par_rv.Cv_iso     = CommandLineArgs::read("-pwcviso_rv", par_rv.Cv_iso,
      "penalidade de taxa das fases isovolumetricas do VD",
      "pvloop por fases");
  par_rv.Cp_fill    = CommandLineArgs::read("-pwcpfill_rv", par_rv.Cp_fill,
      "penalidade de volume do enchimento passivo do VD",
      "pvloop por fases");
  par_rv.Cv_fill    = CommandLineArgs::read("-pwcvfill_rv", par_rv.Cv_fill,
      "penalidade de taxa do enchimento passivo do VD", "pvloop por fases");

  const double cstab = CommandLineArgs::read("-pwcstab", par_lv.C_stab,
      "complacencia de estabilizacao da ejecao, mL/Pa (<=0 = C_cav medida)",
      "pvloop por fases");
  par_lv.C_stab = cstab;
  par_rv.C_stab = cstab;

  const double damp = CommandLineArgs::read("-pwdamp", par_lv.damping,
      "amortecimento do modo automatico das penalidades (0 = deadbeat)",
      "pvloop por fases");
  par_lv.damping = damp;
  par_rv.damping = damp;

  const bool ff = (CommandLineArgs::read("-pwff", par_lv.feedforward ? 1 : 0,
      "1 = antecipacao nas fases isovolumetricas (2/4), compensa o dV/dt da Ta",
      "pvloop por fases") != 0);
  par_lv.feedforward = ff;
  par_rv.feedforward = ff;

  par_lv.fill_rate = CommandLineArgs::read("-pwfillrate", par_lv.fill_rate,
      "vazao de enchimento da fase 5 do VE, mL/ms (<=0 = volta a V_DS)",
      "pvloop por fases");
  par_rv.fill_rate = CommandLineArgs::read("-pwfillrate_rv", par_rv.fill_rate,
      "vazao de enchimento da fase 5 do VD, mL/ms (<=0 = volta a V_DS)",
      "pvloop por fases");

  const double dpmax = CommandLineArgs::read("-pwdpmax", par_lv.dP_max,
      "teto de |dP| por passo, Pa (0 desativa)", "pvloop por fases");
  par_lv.dP_max = dpmax;
  par_rv.dP_max = dpmax;

  // ---- 3. validacoes baratas ------------------------------------------
  if (tim.t_EDP >= tim.RR)
    std::cout << " *** ATENCAO: t_EDP (" << tim.t_EDP << " ms) >= RR ("
              << tim.RR << " ms): a fase 1 ocupa o ciclo inteiro." << std::endl;

  if (par_lv.EDP <= par_lv.P_DS)
    std::cout << " *** ATENCAO: EDP_LV (" << par_lv.EDP << ") <= P_DS_LV ("
              << par_lv.P_DS << "): a contracao atrial ESVAZIA o ventriculo."
              << std::endl;

  if (par_lv.P_ejection <= par_lv.EDP)
    std::cout << " *** ATENCAO: o limiar de ejecao do VE (" << par_lv.P_ejection
              << " Pa) nao supera a EDP (" << par_lv.EDP << " Pa): a fase 2"
              << " termina no primeiro passo e nao havera contracao"
              << " isovolumetrica." << std::endl;

  if (par_rv.R_art > par_lv.R_art)
    std::cout << " *** ATENCAO: R_RV (" << par_rv.R_art << ") > R_LV ("
              << par_lv.R_art << " Pa*ms/mL). A pos-carga pulmonar ficaria"
              << " maior que a sistemica; ver a nota sobre a Tabela 1 em"
              << " piecewise_pvloop.cpp." << std::endl;
}

// ======================================================================
//  Ciclo de vida
// ======================================================================

void PiecewisePVLoop::init(double V_lv0, double V_rv0)
{
  ChamberState * ss[2] = { &st_lv, &st_rv };
  double v0[2] = { V_lv0, V_rv0 };
  ChamberParams * pp[2] = { &par_lv, &par_rv };

  for (int c = 0; c < 2; c++)
  {
    ChamberState & s = *ss[c];
    s.phase     = INIT;
    s.P         = 0.0;
    s.P_art     = pp[c]->P_ejection;
    s.P_prev    = 0.0;
    s.V         = v0[c];
    s.V_prev    = v0[c];
    s.Vdot      = 0.0;
    s.V_lock    = v0[c];
    s.V_DS      = (pp[c]->DSV > 0.0) ? pp[c]->DSV : v0[c];
    s.V_ED      = v0[c];
    s.V_ES      = v0[c];
    s.V_ref     = v0[c];
    s.has_ejected = false;
    s.t_phase   = 0.0;
    s.n_clipped = 0;

    // Ainda nao ha nada medido. Um valor positivo de partida evita divisao
    // por zero se o modo automatico for consultado antes do fim da fase 0;
    // ele e substituido pela medida real assim que ela existir.
    s.C_cav     = (v0[c] > 0.0) ? 0.25e-3 * v0[c] : 1.0e-2;   // [mL/Pa]
  }
}

int PiecewisePVLoop::beat_of(double t) const
{
  if (t < tim.t_init) return -1;
  return static_cast<int>(std::floor((t - tim.t_init) / tim.RR));
}

double PiecewisePVLoop::clip_dP(double dP, ChamberState & st,
                                const ChamberParams & p) const
{
  if (p.dP_max <= 0.0) return dP;
  if (dP >  p.dP_max) { st.n_clipped++; return  p.dP_max; }
  if (dP < -p.dP_max) { st.n_clipped++; return -p.dP_max; }
  return dP;
}

double PiecewisePVLoop::compliance_from(double V_hi, double P_hi,
                                        double V_lo, double P_lo)
{
  const double dV = V_hi - V_lo;
  const double dP = P_hi - P_lo;

  // Uma cavidade que ganhou volume ao ser despressurizada, ou que nao
  // respondeu, nao produz uma complacencia utilizavel. Devolver 0 sinaliza
  // "sem medida" e deixa o chamador manter a estimativa anterior, o que e
  // muito melhor do que propagar um ganho negativo para o controlador.
  if (dP <= 1e-6 || dV <= 0.0) return 0.0;   // dP em Pa
  return dV / dP;
}

double PiecewisePVLoop::penalty_dP(ChamberState & st, const ChamberParams & p,
                                   double V_target, double Cp, double Cv,
                                   double dt, bool allow_ff) const
{
  double dP;

  if (Cp > 0.0)
  {
    // ---- formula literal do artigo ----------------------------------
    // dP = -(1/Cp)(V - V_target) - (1/Cv) dV/dt, com o fator g = dt no
    // modo de taxa. Cv <= 0 desliga so o termo derivativo.
    const double g = rate_form ? dt : 1.0;
    dP = -(st.V - V_target) / Cp;
    if (Cv > 0.0) dP -= st.Vdot / Cv;
    dP *= g;
  }
  else
  {
    // ---- modo automatico --------------------------------------------
    // dP = -[ (V - V_target) + damping * dt * dV/dt ] / C_cav.
    // Cp = C_cav e o ganho deadbeat: um passo corrigiria o desvio inteiro
    // se a cavidade fosse linear com a complacencia medida. Ver a nota no
    // cabecalho de piecewise_pvloop.hpp.
    const double C = (st.C_cav > 1e-12) ? st.C_cav : 1e-12;
    if (allow_ff && p.feedforward && dt > 0.0)
    {
      // Antecipacao: a parte de dV/dt que nao veio da pressao do passo
      // anterior (a Ta) deve continuar no proximo passo; compensa ja.
      const double Vdot_act = st.Vdot - C * (st.P - st.P_prev) / dt;
      dP = -((st.V - V_target) + dt * Vdot_act) / C;
    }
    else
      dP = -((st.V - V_target) + p.damping * dt * st.Vdot) / C;
  }

  return clip_dP(dP, st, p);
}

void PiecewisePVLoop::arterial_runoff(ChamberState & st,
                                      const ChamberParams & p,
                                      double dt) const
{
  // Escoamento diastolico da arteria: C dP/dt = -P/R, sem fluxo de entrada.
  // So faz sentido no modo continuo; no modo literal P_art e reiniciado em
  // P_ejection a cada abertura da valva e nao tem historia entre batimentos.
  if (!continuous_wk || p.C_art <= 0.0 || p.R_art <= 0.0) return;
  if (wk_hold_until_ejection && !st.has_ejected) return;
  st.P_art += dt * (-st.P_art / p.R_art) / p.C_art;
  if (st.P_art < 0.0) st.P_art = 0.0;
}

void PiecewisePVLoop::advance_one(ChamberState & st, const ChamberParams & p,
                                  double t_end, double dt) const
{
  // O fator dt/1 da ambiguidade literal-vs-taxa vive dentro de
  // penalty_dP(), que e o unico lugar onde ele tem efeito.
  const double P_old = st.P;

  switch (st.phase)
  {
    // ---- Fase 0: rampa de inflacao ate a diastase --------------------
    case INIT:
    {
      const double s = (tim.t_init > 0.0) ? clamp01(t_end / tim.t_init) : 1.0;
      st.P = p.P_DS * s;
      arterial_runoff(st, p, dt);
      break;
    }

    // ---- Fase 1: enchimento diastolico ativo -------------------------
    case FILLING:
    {
      // A rampa e ancorada no INSTANTE DE ENTRADA na fase, nao no inicio
      // nominal do batimento. Assim ela continua correta quando a guarda de
      // relogio force a entrada com alguns passos de atraso, e a rampa e a
      // condicao de saida abaixo usam a mesma referencia.
      const double s = (tim.t_EDP > 0.0)
                     ? clamp01((t_end - st.t_phase) / tim.t_EDP) : 1.0;
      st.P = p.P_DS + (p.EDP - p.P_DS) * s;
      arterial_runoff(st, p, dt);
      break;
    }

    // ---- Fases 2 e 4: isovolumetricas --------------------------------
    //
    // Controlador proporcional-derivativo sobre a pressao, com o volume
    // travado em V_lock. O termo proporcional corrige o desvio acumulado;
    // o derivativo amortece. Sem o segundo, a pressao oscila com periodo 2
    // -- e um esquema explicito acoplado a uma cavidade muito rigida.
    case IVC:
    case IVR:
    {
      st.P += penalty_dP(st, p, st.V_lock, p.Cp_iso, p.Cv_iso, dt);
      arterial_runoff(st, p, dt);
      break;
    }

    // ---- Fase 3: ejecao, Windkessel de tres elementos ----------------
    //
    //     C_art dP_art/dt + P_art/R_art = Q,   P = P_art + R_c Q
    //
    // com Q = -dV/dt o fluxo que SAI do ventriculo, positivo durante a
    // ejecao. Sempre integrada como equacao diferencial, nas duas leituras:
    // aqui o artigo escreve a derivada explicitamente. R_c = 0 e o Windkessel
    // de dois elementos do artigo (P = P_art).
    case EJECTION:
    {
      if (p.C_art > 0.0 && p.R_art > 0.0 && dt > 0.0)
      {
        // Complacencia de estabilizacao. Ver a nota da fase 3 no
        // cabecalho da classe: ela some no estado estacionario.
        const double C_s = (p.C_stab > 0.0) ? p.C_stab
                         : ((st.C_cav > 0.0) ? st.C_cav : 0.0);

        // Pressao REAL da cavidade no fim do passo anterior e sua variacao.
        // No primeiro passo apos a abertura e a pressao da fase 2, que ja
        // passou do limiar e vinha subindo: usar P_art aqui faria a pressao
        // cair de uma vez e o ventriculo despejar o volume num unico passo.
        const double P_n     = st.P;
        const double dP_prev = st.P - st.P_prev;

        // Parte de dV/dt que NAO vem da variacao de pressao do proprio
        // passo anterior, isto e, o encurtamento imposto pela tensao
        // ativa. E ela que alimenta a arteria.
        const double Vdot_act = st.Vdot - C_s * dP_prev / dt;

        // Euler implicito, resolvido junto para P_art+ e Q:
        //   C_art (P_art+ - P_art)/dt + P_art+/R = Q
        //   Q = -Vdot_act - C_s (P+ - P_n)/dt,   P+ = P_art+ + R_c Q
        // A primeira da P_art+ = a + b Q; substituindo na segunda sai Q.
        // Com R_c = 0 reduz a  C_tot (P+ - P)/dt + P+/R = -Vdot_act.
        const double b = 1.0 / (p.C_art / dt + 1.0 / p.R_art);
        const double a = b * p.C_art * st.P_art / dt;
        const double k = C_s / dt;
        const double Q     = (-Vdot_act - k * (a - P_n)) / (1.0 + k * (b + p.R_c));
        const double P_new = a + b * Q;

        st.P_art     += clip_dP(P_new - st.P_art, st, p);
        if (st.P_art < 0.0) st.P_art = 0.0;

        // Valva arterial aberta: a cavidade ve a arteria atraves de R_c.
        st.P = st.P_art + p.R_c * Q;
      }
      else
        st.P = st.P_art;
      break;
    }

    // ---- Fase 5: enchimento passivo ----------------------------------
    case PASSIVE:
    {
      // Vazao prescrita: rampa a partir do volume de entrada (V_lock).
      const double V_target = (p.fill_rate > 0.0)
                            ? st.V_lock + p.fill_rate * (t_end - st.t_phase)
                            : st.V_DS;
      // Sem antecipacao: aqui a Ta ja e nula e, em pressao baixa, a cavidade
      // e complacente demais para ela (oscila com periodo 2).
      st.P += penalty_dP(st, p, V_target, p.Cp_fill, p.Cv_fill, dt, false);
      arterial_runoff(st, p, dt);
      break;
    }
  }

  st.P_prev = P_old;
}

void PiecewisePVLoop::advance_pressure(double t_end, double dt)
{
  advance_one(st_lv, par_lv, t_end, dt);
  advance_one(st_rv, par_rv, t_end, dt);
}

// ======================================================================
//  Transicoes
// ======================================================================

void PiecewisePVLoop::transition_one(ChamberState & st, const ChamberParams & p,
                                     double t, const char * tag)
{
  const double in_phase = t - st.t_phase;

  // ---- guarda de relogio ---------------------------------------------
  //
  // O inicio de cada batimento reinicia a fase 1, aconteca o que
  // acontecer. Sem isso, uma cavidade que nunca alcanca o limiar de ejecao
  // (Ta baixa demais) ou que nunca ve o fluxo reverter fica presa numa
  // fase para sempre e a rodada inteira sai sem sentido, silenciosamente.
  const int b = beat_of(t);
  if (b >= 0)
  {
    const double t0 = filling_start(b);
    const double in_filling_window = t - t0;

    if (in_filling_window >= 0.0 && in_filling_window < tim.t_EDP
        && st.phase != FILLING)
    {
      if (st.phase != PASSIVE && st.phase != INIT)
        std::cout << "  [PW " << tag << "] *** batimento " << b + 1
                  << " comecou com a cavidade ainda em "
                  << phase_name(st.phase)
                  << "; forcada para FILLING pelo relogio." << std::endl;

      // O volume de diastase automatico e fixado no fim da fase 0.
      if (st.phase == INIT)
      {
        if (p.DSV <= 0.0) st.V_DS = st.V;

        // Primeira estimativa da complacencia: a malha foi de (V_ref, 0)
        // ate (V, P_DS) durante a rampa de inflacao. Grosseira, porque a
        // curva passiva e exponencial e esta secante e longa demais, mas
        // e a unica disponivel antes do primeiro fim de diastole e ja
        // poe o ganho na ordem de grandeza certa.
        const double c = compliance_from(st.V, p.P_DS, st.V_ref, 0.0);
        if (c > 0.0) st.C_cav = c;

        std::cout << "  [PW " << tag << "] fim da fase 0: V_diastase = "
                  << st.V_DS << " mL, C_cav ~ " << st.C_cav
                  << " mL/Pa" << std::endl;
      }

      st.phase   = FILLING;
      st.t_phase = t;
      return;
    }
  }

  switch (st.phase)
  {
    case INIT:
      // A saida da fase 0 e tratada pela guarda de relogio acima, que
      // dispara assim que t alcanca t_init.
      break;

    case FILLING:
      // Fim da diastole. Aqui se registram V_ED e o volume travado da
      // contracao isovolumetrica, e e aqui que a ativacao eletrica ja
      // deve ter sido disparada (ver activation_time()).
      if (in_phase >= tim.t_EDP - 1e-9)
      {
        st.V_ED   = st.V;
        st.V_lock = st.V;

        // Secante da curva passiva entre diastase e fim de diastole. E a
        // medida boa: os dois pontos sao proximos, estao ambos na perna
        // diastolica e foram produzidos por esta malha com este material.
        // E ela que passa a alimentar o ganho automatico das fases 2/4/5.
        const double c = compliance_from(st.V_ED, p.EDP, st.V_DS, p.P_DS);
        if (c > 0.0) st.C_cav = c;

        st.phase  = IVC;
        st.t_phase = t;
        std::cout << "  [PW " << tag << "] fim da diastole: V_ED = "
                  << st.V_ED << " mL, P = " << st.P << " Pa, C_cav = "
                  << st.C_cav << " mL/Pa -> IVC" << std::endl;
      }
      break;

    case IVC:
    {
      // Abertura da valva arterial.
      const double thr = continuous_wk ? st.P_art : p.P_ejection;
      if (st.P >= thr)
      {
        if (!continuous_wk) st.P_art = p.P_ejection;
        st.phase   = EJECTION;
        st.has_ejected = true;
        st.t_phase = t;
        std::cout << "  [PW " << tag << "] abertura da valva: P = " << st.P
                  << " Pa >= " << thr << " Pa -> EJECTION" << std::endl;
      }
      break;
    }

    case EJECTION:
    {
      // Criterio do artigo: reversao do fluxo. O guarda de duracao minima
      // existe porque, no primeiro passo apos a abertura, dV/dt ainda e o
      // valor herdado da fase isovolumetrica -- proximo de zero e com
      // sinal decidido pelo ruido do solver.
      const bool reversed = (st.Vdot > 0.0);

      // Rede de seguranca. Se a arteria drena mais rapido do que o
      // ventriculo relaxa, dV/dt nunca fica positivo e a cavidade fica
      // presa na ejecao ate o relogio do proximo batimento -- sem V_ES,
      // sem fases 4 e 5, e sem nenhum aviso. Isso acontece quando R_art e
      // baixa demais ou C_art pequena demais para a complacencia desta
      // malha, entao vale a pena dizer isso em voz alta.
      const bool collapsed = (st.P <= p.P_atrium);

      if (in_phase >= tim.min_phase && (reversed || collapsed))
      {
        if (!reversed)
          std::cout << "  [PW " << tag << "] *** o fluxo NAO reverteu; a"
                    << " ejecao foi encerrada porque P caiu a " << st.P
                    << " Pa <= P_atrio (" << p.P_atrium << "). Aumente"
                    << " R_art/C_art ou reveja o relaxamento." << std::endl;

        st.V_ES   = st.V;
        st.V_lock = st.V;
        st.phase  = IVR;
        st.t_phase = t;
        std::cout << "  [PW " << tag << "] fim da ejecao: V_ES = " << st.V_ES
                  << " mL, P = " << st.P << " Pa -> IVR" << std::endl;
      }
      break;
    }

    case IVR:
      // Abertura da valva atrioventricular.
      if (in_phase >= tim.min_phase && st.P <= p.P_atrium)
      {
        st.phase   = PASSIVE;
        st.t_phase = t;
        st.V_lock  = st.V;          // origem da rampa de enchimento
        std::cout << "  [PW " << tag << "] P = " << st.P << " Pa <= P_atrio ("
                  << p.P_atrium << ") -> PASSIVE" << std::endl;
      }
      break;

    case PASSIVE:
      // Sai apenas pelo relogio, na guarda no topo desta funcao.
      break;
  }
}

void PiecewisePVLoop::update_volumes(double t, double dt,
                                     double V_lv, double V_rv)
{
  ChamberState * ss[2] = { &st_lv, &st_rv };
  ChamberParams * pp[2] = { &par_lv, &par_rv };
  double vv[2] = { V_lv, V_rv };
  const char * tag[2] = { "LV", "RV" };

  for (int c = 0; c < 2; c++)
  {
    ChamberState & s = *ss[c];
    s.V_prev = s.V;
    s.V      = vv[c];
    s.Vdot   = (dt > 0.0) ? (s.V - s.V_prev) / dt : 0.0;

    transition_one(s, *pp[c], t, tag[c]);
  }
}

// ======================================================================
//  Consultas e saida
// ======================================================================

double PiecewisePVLoop::ejection_fraction_LV() const
{
  if (st_lv.V_ED <= 0.0) return 0.0;
  return 100.0 * (st_lv.V_ED - st_lv.V_ES) / st_lv.V_ED;
}

double PiecewisePVLoop::ejection_fraction_RV() const
{
  if (st_rv.V_ED <= 0.0) return 0.0;
  return 100.0 * (st_rv.V_ED - st_rv.V_ES) / st_rv.V_ED;
}

std::string PiecewisePVLoop::describe() const
{
  std::ostringstream os;

  os << " PV loop por fases (Wang et al. 2025, eLife 14:RP106555)\n";
  os << "   unidades: mL, Pa, ms\n";
  os << "   periodo RR = " << tim.RR << " ms ("
     << (tim.RR > 0.0 ? 60000.0 / tim.RR : 0.0) << " bpm)\n";
  os << "   fase 0 (inflacao) = " << tim.t_init << " ms;"
     << " fase 1 (chute atrial) = " << tim.t_EDP << " ms\n";
  os << "   ativacao eletrica do primeiro batimento em t = "
     << activation_time(0) << " ms\n";
  os << "   fases 2/4/5: " << (rate_form ? "forma de TAXA (P += dt*dP)"
                                         : "forma LITERAL do artigo (P += dP)")
     << "\n";
  os << "   valva arterial: " << (!continuous_wk
        ? "limiar fixo P_ejection, P_art reiniciado a cada ejecao"
        : wk_hold_until_ejection
        ? "Windkessel continuo, limiar = P_art corrente, drena so apos a 1a ejecao"
        : "Windkessel continuo, limiar = P_art corrente") << "\n";

  const ChamberParams * pp[2] = { &par_lv, &par_rv };
  const char * tag[2] = { "VE", "VD" };

  for (int c = 0; c < 2; c++)
  {
    const ChamberParams & p = *pp[c];
    os << "   " << tag[c] << ": P_DS = " << p.P_DS
       << ", EDP = " << p.EDP
       << ", P_ejecao = " << p.P_ejection
       << ", P_atrio = " << p.P_atrium << " Pa\n";
    os << "       C_art = " << p.C_art << " mL/Pa, R_art = " << p.R_art
       << " Pa*ms/mL, R_c = " << p.R_c << " Pa*ms/mL, DSV = "
       << (p.DSV > 0.0 ? std::to_string(p.DSV) : std::string("auto")) << "\n";
    os << "       penalidade isovolumetrica: "
       << (p.Cp_iso > 0.0
             ? "Cp_iso = " + std::to_string(p.Cp_iso) + ", Cv_iso = "
               + std::to_string(p.Cv_iso)
             : "AUTOMATICA (C_cav medido, damping = "
               + std::to_string(p.damping) + ")") << "\n";
    os << "       penalidade de enchimento:  "
       << (p.Cp_fill > 0.0
             ? "Cp_fill = " + std::to_string(p.Cp_fill) + ", Cv_fill = "
               + std::to_string(p.Cv_fill)
             : "AUTOMATICA (C_cav medido, damping = "
               + std::to_string(p.damping) + ")") << "\n";
    os << "       fase 5: "
       << (p.fill_rate > 0.0 ? "vazao prescrita " + std::to_string(p.fill_rate) + " mL/ms"
                             : std::string("PD de volta a V_DS"))
       << (p.feedforward ? "; fases 2/4 com antecipacao (-pwff)" : "") << "\n";
  }

  os << "   NOTA: as penalidades Cp/Cv NAO sao tabeladas no artigo"
     << " (aparecem apenas como \"inverse of penalty terms\").\n";
  os << "         No modo automatico o ganho e a complacencia diastolica"
     << " medida na propria malha,\n";
  os << "         o que garante estabilidade; um Cp fixo so converge se"
     << " C_cav/Cp < 2.";

  return os.str();
}

std::string PiecewisePVLoop::history_header()
{
  return "# t[ms] phase_LV P_LV[Pa] V_LV[mL] dVdt_LV[mL/ms] Part_LV[Pa]"
         " phase_RV P_RV[Pa] V_RV[mL] dVdt_RV[mL/ms] Part_RV[Pa]";
}

std::string PiecewisePVLoop::history_row(double t) const
{
  std::ostringstream os;
  os << std::scientific << std::setprecision(8);
  os << t
     << " " << static_cast<int>(st_lv.phase)
     << " " << st_lv.P << " " << st_lv.V << " " << st_lv.Vdot
     << " " << st_lv.P_art
     << " " << static_cast<int>(st_rv.phase)
     << " " << st_rv.P << " " << st_rv.V << " " << st_rv.Vdot
     << " " << st_rv.P_art;
  return os.str();
}