#include "cell_warmup.hpp"
#include "cellmodel.hpp"
#include "util/command_line_args.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

// ===========================================================================
//  Opcoes
// ===========================================================================

CellWarmup::Options
CellWarmup::options_from_command_line(double default_bcl_ms,
                                      double default_dt_ms,
                                      const std::string & model_name,
                                      const std::string & basename)
{
  Options o;

  o.enabled = (CommandLineArgs::read("-warmup", 0) != 0);
  if (!o.enabled) return o;

  o.beats       = CommandLineArgs::read("-warmup_beats", 100);
  o.bcl_ms      = CommandLineArgs::read("-warmup_bcl", default_bcl_ms);
  o.dt_ms       = CommandLineArgs::read("-warmup_dt", default_dt_ms);
  o.stim_amp    = CommandLineArgs::read("-warmup_stimamp",
                                        CommandLineArgs::read("-stimamp", -53.0));
  o.stim_dur_ms = CommandLineArgs::read("-warmup_stimdur", 2.0);
  o.phase_ms    = CommandLineArgs::read("-warmup_phase", 0.0);
  o.phase_from_lat = (CommandLineArgs::read("-warmup_lat", 0) != 0);
  o.sample_ms   = CommandLineArgs::read("-warmup_sample", 1.0);
  o.ab_bins     = CommandLineArgs::read("-warmup_abbins", 1);
  o.lambda      = CommandLineArgs::read("-warmup_lambda", 1.0);
  o.tol         = CommandLineArgs::read("-warmup_tol", 0.0);
  o.model_name  = model_name;

  // Nome default do cache derivado da malha, para que rodadas diferentes na
  // mesma malha reaproveitem o mesmo arquivo. "none" desliga o cache.
  std::string base = basename;
  const std::size_t pos = base.find(".xml");
  if (pos != std::string::npos) base = base.substr(0, pos);
  o.cache = CommandLineArgs::read("-warmup_cache", base + "_warmup.dat");
  if (o.cache == "none" || o.cache == "off") o.cache.clear();

  return o;
}

// ===========================================================================
//  Construcao
// ===========================================================================

CellWarmup::CellWarmup(CellModel * m, const Options & o)
  : model(m), opt(o)
{
  nvars = model->get_num_state_vars();
  unit  = model->native_time_unit_ms();
  if (unit <= 0.0) unit = 1.0;

  // Fase negativa conta de tras para frente: -100 significa "100 ms antes do
  // repouso", isto e, 100 ms antes do proximo estimulo.
  phase_wrapped = std::fmod(opt.phase_ms, opt.bcl_ms);
  if (phase_wrapped < 0.0) phase_wrapped += opt.bcl_ms;
}

const char * CellWarmup::type_name(int t)
{
  switch (t)
  {
    case EPI:   return "epi";
    case MCELL: return "mid";
    case ENDO:  return "endo";
    case APEX:  return "apex";
    case BASE:  return "base";
    default:    return "?";
  }
}

int CellWarmup::group_of(int cell) const
{
  if (cell < 0 || cell >= (int) cell_group.size()) return 0;
  return cell_group[cell];
}

// ===========================================================================
//  Agrupamento
//
//  O ciclo limite depende so dos parametros da celula, nao da sua posicao.
//  Duas celulas com o mesmo tipo e a mesma coordenada apicobasal tem
//  exatamente o mesmo estado estacionario, entao basta pre-condicionar uma
//  celula por combinacao (tipo, faixa de ab).
// ===========================================================================

void CellWarmup::build_groups(const arma::ivec & types, const arma::vec & ab,
                              int n_cells)
{
  // Valor NEUTRO da coordenada apicobasal: no ToRORd-Land o gradiente e
  // GKs * 0.2^(2*ab - 1), que vale exatamente 1 em ab = 0.5. Usar esse valor
  // quando o gradiente e ignorado deixa o warm up igual ao ciclo limite do
  // modelo original, e -- por nao depender da malha -- torna o cache
  // reaproveitavel entre malhas diferentes.
  const double AB_NEUTRAL = 0.5;

  const bool has_types = (types.n_elem == (arma::uword) n_cells);
  const int  nbins     = std::max(1, opt.ab_bins);

  // Com uma unica faixa o gradiente apicobasal simplesmente NAO entra no
  // agrupamento: sobra um grupo por tipo celular (tipicamente 3). O erro
  // que isso introduz e da ordem de 1% nas variaveis lentas -- duas ordens
  // de grandeza menor que o erro das condicoes iniciais publicadas, que e o
  // que o warm up existe para corrigir. Use -warmup_abbins 8 se quiser
  // resolver o gradiente.
  const bool use_ab = (ab.n_elem == (arma::uword) n_cells) && (nbins > 1);

  cell_group.assign(n_cells, 0);
  grp_type.clear();
  grp_ab.clear();

  std::map<long, int>  seen;
  std::vector<double>  ab_sum;
  std::vector<int>     ab_cnt;

  for (int i = 0; i < n_cells; i++)
  {
    const int t = has_types ? (int) types(i) : model->get_celltype();

    int b = 0;
    double a = AB_NEUTRAL;
    if (use_ab)
    {
      a = ab(i);
      double ac = a;
      if (ac < 0.0) ac = 0.0;
      else if (ac > 1.0) ac = 1.0;
      b = (int) (ac * nbins);
      if (b >= nbins) b = nbins - 1;
    }

    const long key = (long) t * 1000L + (long) b;

    std::map<long,int>::iterator it = seen.find(key);
    int g;
    if (it == seen.end())
    {
      g = (int) grp_type.size();
      seen[key] = g;
      grp_type.push_back(t);
      ab_sum.push_back(0.0);
      ab_cnt.push_back(0);
    }
    else
      g = it->second;

    cell_group[i] = g;
    ab_sum[g] += a;
    ab_cnt[g] += 1;
  }

  // A ab representativa do grupo e a MEDIA das celulas que caem nele, nao o
  // centro da faixa: com poucas faixas isso reduz bastante o erro cometido
  // ao tratar um gradiente continuo como discreto.
  grp_ab.resize(grp_type.size());
  for (std::size_t g = 0; g < grp_type.size(); g++)
    grp_ab[g] = (ab_cnt[g] > 0) ? ab_sum[g] / ab_cnt[g] : 0.5;
}

// ===========================================================================
//  Integracao de um trecho do protocolo
// ===========================================================================

void CellWarmup::pace(double * y, double length_ms, arma::mat * rec,
                      bool fill_times)
{
  if (length_ms <= 0.0)
  {
    // Fase 0 sem amostragem: o estado ja e o pedido.
    if (rec)
    {
      rec->set_size(nvars, 1);
      rec->col(0) = arma::vec(y, nvars);
      if (fill_times) { snap_t.set_size(1); snap_t(0) = 0.0; }
    }
    return;
  }

  // O passo efetivo divide exatamente o intervalo, para que o ultimo passo
  // caia sobre length_ms e a fase pedida seja atingida sem sobra.
  int nsteps = (int) std::llround(length_ms / opt.dt_ms);
  if (nsteps < 1) nsteps = 1;
  const double dt_ms_eff = length_ms / nsteps;
  const double dt_n      = dt_ms_eff / unit;

  int nsub = 1;
  if (rec)
  {
    nsub = (int) std::llround(opt.sample_ms / dt_ms_eff);
    if (nsub < 1) nsub = 1;
  }

  std::vector<arma::vec> snaps;
  std::vector<double>    times;

  for (int k = 0; k < nsteps; k++)
  {
    const double tk_ms = k * dt_ms_eff;

    if (rec && (k % nsub == 0))
    {
      snaps.push_back(arma::vec(y, nvars));
      times.push_back(tk_ms);
    }

    // Estimulo no INICIO do intervalo. tk_ms e o instante em que o passo
    // comeca, entao a comparacao e com o inicio, nao com o fim, do passo.
    const double istim = (tk_ms < opt.stim_dur_ms) ? opt.stim_amp : 0.0;

    model->advance(y, tk_ms / unit, dt_n, istim);
  }

  if (rec)
  {
    // A ultima amostra e sempre o fim do intervalo, mesmo que nao caia num
    // multiplo de sample_ms: e ela que fecha o ciclo em [0, BCL].
    snaps.push_back(arma::vec(y, nvars));
    times.push_back(length_ms);

    rec->set_size(nvars, snaps.size());
    for (std::size_t j = 0; j < snaps.size(); j++)
      rec->col(j) = snaps[j];

    if (fill_times)
    {
      snap_t.set_size(times.size());
      for (std::size_t j = 0; j < times.size(); j++) snap_t(j) = times[j];
    }
  }
}

double CellWarmup::drift(const arma::vec & a, const arma::vec & b)
{
  double m = 0.0;
  for (arma::uword i = 0; i < a.n_elem; i++)
  {
    const double den = std::max(std::fabs(a(i)), std::fabs(b(i)));
    // Variaveis numericamente nulas (Jrelnp ~ 1e-24 em repouso) nao dizem
    // nada sobre convergencia e dominariam qualquer medida relativa.
    if (den < 1.0e-8) continue;
    m = std::max(m, std::fabs(a(i) - b(i)) / den);
  }
  return m;
}

// ===========================================================================
//  Execucao
// ===========================================================================

bool CellWarmup::run(const arma::ivec & types, const arma::vec & ab,
                     int n_cells)
{
  if (!opt.enabled || model == nullptr || n_cells <= 0) return false;

  if (opt.bcl_ms <= 0.0 || opt.dt_ms <= 0.0)
  {
    std::cout << " *** warm up: invalid BCL (" << opt.bcl_ms << " ms) or dt ("
              << opt.dt_ms << " ms); pre-conditioning"
              << " disabled." << std::endl;
    return false;
  }

  std::cout << "\n=== 0D cell model warm up ===" << std::endl;

  build_groups(types, ab, n_cells);

  const int ngroups = (int) grp_type.size();
  const long steps_per_beat = std::llround(opt.bcl_ms / opt.dt_ms);

  std::cout << " Model: " << (opt.model_name.empty() ? "?" : opt.model_name)
            << " (" << nvars << " state variables, native unit "
            << unit << " ms)" << std::endl;
  std::cout << " Protocol: " << opt.beats << " beats, BCL = "
            << opt.bcl_ms << " ms, dt = " << opt.dt_ms << " ms" << std::endl;
  std::cout << " Stimulus: amplitude " << opt.stim_amp << ", duration "
            << opt.stim_dur_ms << " ms" << std::endl;
  std::cout << " Phase of the stored state: " << phase_wrapped
            << " ms after the stimulus";
  if (phase_wrapped == 0.0) std::cout << " (rest)";
  std::cout << std::endl;
  if (opt.phase_from_lat)
    std::cout << "   phase shifted per node: phase_i = " << phase_wrapped
              << " - lat_i (sampled every " << opt.sample_ms << " ms)"
              << std::endl;
  std::cout << " Stretch during warm up: lambda = " << opt.lambda
            << ", d(lambda)/dt = 0" << std::endl;
  std::cout << " Cells: " << n_cells << " node(s) -> " << ngroups
            << " group(s) (type x ab bin, -warmup_abbins "
            << opt.ab_bins << ")" << std::endl;
  if (opt.ab_bins <= 1 && ab.n_elem == (arma::uword) n_cells)
    std::cout << "   apicobasal gradient IGNORED during warm up: one group per"
              << " cell type, ab = 0.5 (neutral). The gradient stays"
              << " active in the simulation." << std::endl;
  std::cout << " Cost: " << ngroups << " x " << (opt.beats + 1) << " x "
            << steps_per_beat << " = "
            << (double) ngroups * (opt.beats + 1) * steps_per_beat
            << " ODE steps" << std::endl;

  if (model->lat_var_index() >= 0)
    std::cout << " *** WARNING: this model is phenomenological and is triggered"
              << " by activation time, not by a stimulus current."
              << " The warm up probably makes no sense for it."
              << std::endl;

  // ---- cache --------------------------------------------------------
  if (load_cache())
  {
    is_ready = true;
    std::cout << " State read from cache " << opt.cache
              << "; pre-conditioning skipped." << std::endl;
    std::cout << "=== Warm up done ===\n" << std::endl;
    return true;
  }

  // ---- model state to preserve --------------------------------------
  const int    type0   = model->get_celltype();
  const double ab0     = model->get_apicobasal();
  const double lam0    = model->get_stretch();
  const double lamdot0 = model->get_stretch_rate();

  rest_state.assign(ngroups, arma::vec(nvars, arma::fill::zeros));
  traj.assign(ngroups, arma::mat());

  arma::vec y(nvars), y_prev(nvars);

  for (int g = 0; g < ngroups; g++)
  {
    model->set_celltype(grp_type[g]);
    model->set_apicobasal(grp_ab[g]);
    model->set_stretch(opt.lambda);
    model->set_stretch_rate(0.0);

    // Condicoes iniciais publicadas do modelo, ja para o tipo do grupo.
    model->init(y.memptr());

    std::cout << " [grupo " << g + 1 << "/" << ngroups << "] "
              << type_name(grp_type[g]) << ", ab = "
              << std::fixed << std::setprecision(3) << grp_ab[g]
              << std::defaultfloat << " ... " << std::flush;

    double d = 0.0;
    int    beats_run = 0;

    for (int beat = 0; beat < opt.beats; beat++)
    {
      y_prev = y;
      pace(y.memptr(), opt.bcl_ms, nullptr, false);
      beats_run++;

      if (!y.is_finite())
      {
        std::cout << "\n *** warm up: the model diverged (NaN/Inf) on"
                  << " beat " << beat + 1 << ". Reduce -warmup_dt."
                  << std::endl;
        model->set_celltype(type0);
        model->set_apicobasal(ab0);
        model->set_stretch(lam0);
        model->set_stretch_rate(lamdot0);
        return false;
      }

      d = drift(y, y_prev);

      // Parada antecipada: o ciclo limite ja repete dentro da tolerancia.
      if (opt.tol > 0.0 && d < opt.tol) break;
    }

    // ---- ultimo trecho, ate a fase pedida --------------------------
    if (opt.phase_from_lat)
    {
      // Precisa do batimento inteiro para poder amostrar qualquer fase.
      pace(y.memptr(), opt.bcl_ms, &traj[g], (g == 0));
    }
    else if (phase_wrapped > 0.0)
    {
      // So o pedaco inicial do batimento seguinte: assim a fase e atingida
      // exatamente, sem interpolacao.
      pace(y.memptr(), phase_wrapped, nullptr, false);
    }

    rest_state[g] = y;

    std::cout << beats_run << " beat(s), final drift "
              << std::scientific << std::setprecision(2) << d
              << std::defaultfloat;
    if (nvars > 0)
      std::cout << ", V = " << std::fixed << std::setprecision(2) << y(0)
                << " mV" << std::defaultfloat;
    std::cout << std::endl;

    if (opt.tol > 0.0 && d >= opt.tol)
      std::cout << "   *** WARNING: drift still above -warmup_tol ("
                << opt.tol << ") after " << beats_run
                << " beats; increase -warmup_beats." << std::endl;
  }

  // ---- restore the model state --------------------------------------
  model->set_celltype(type0);
  model->set_apicobasal(ab0);
  model->set_stretch(lam0);
  model->set_stretch_rate(lamdot0);

  is_ready = true;
  save_cache();

  std::cout << "=== Warm up done ===\n" << std::endl;
  return true;
}

// ===========================================================================
//  Consulta
// ===========================================================================

void CellWarmup::state_for_cell(int cell, double lat_ms, double * y) const
{
  if (!is_ready || y == nullptr) return;

  const int g = group_of(cell);

  if (!opt.phase_from_lat || traj.empty() || traj[g].n_cols == 0)
  {
    std::memcpy(y, rest_state[g].memptr(), sizeof(double) * nvars);
    return;
  }

  // Fase propria deste no: a fase pedida, atrasada pelo seu tempo de
  // ativacao. Um no que ativa em lat_ms deve estar, no instante 0 da
  // simulacao, lat_ms ANTES do proprio estimulo.
  double ph = std::fmod(phase_wrapped - lat_ms, opt.bcl_ms);
  if (ph < 0.0) ph += opt.bcl_ms;

  const arma::uword n = snap_t.n_elem;
  if (n == 0) { std::memcpy(y, rest_state[g].memptr(), sizeof(double)*nvars); return; }

  // Busca binaria: primeira amostra com tempo > ph.
  const double * t0 = snap_t.memptr();
  const double * up = std::upper_bound(t0, t0 + n, ph);
  arma::uword k = (arma::uword) (up - t0);
  if (k == 0) k = 1;
  if (k >= n) k = n - 1;

  const double ta = snap_t(k - 1), tb = snap_t(k);
  const double w  = (tb > ta) ? (ph - ta) / (tb - ta) : 0.0;

  const arma::mat & M = traj[g];
  for (int i = 0; i < nvars; i++)
    y[i] = (1.0 - w) * M(i, k - 1) + w * M(i, k);
}

// ===========================================================================
//  Cache
//
//  A assinatura contem TODOS os parametros que determinam o resultado,
//  inclusive a lista de grupos. Qualquer diferenca invalida o cache: e
//  preferivel recalcular a arriscar comecar uma simulacao do estado
//  estacionario de outra configuracao.
// ===========================================================================

std::string CellWarmup::fingerprint() const
{
  std::ostringstream s;
  s << std::setprecision(10);
  s << "model=" << opt.model_name
    << " nvars=" << nvars
    << " beats=" << opt.beats
    << " bcl=" << opt.bcl_ms
    << " dt=" << opt.dt_ms
    << " amp=" << opt.stim_amp
    << " sdur=" << opt.stim_dur_ms
    << " phase=" << phase_wrapped
    << " lat=" << (opt.phase_from_lat ? 1 : 0)
    << " sample=" << opt.sample_ms
    << " lam=" << opt.lambda
    << " tol=" << opt.tol
    << " groups=" << grp_type.size();

  for (std::size_t g = 0; g < grp_type.size(); g++)
    s << " [" << grp_type[g] << "," << std::fixed << std::setprecision(6)
      << grp_ab[g] << std::defaultfloat << "]";

  return s.str();
}

bool CellWarmup::load_cache()
{
  if (opt.cache.empty()) return false;

  std::ifstream in(opt.cache.c_str());
  if (!in.is_open()) return false;

  std::string line;
  if (!std::getline(in, line)) return false;
  if (line.rfind("# cardiax cell warm-up cache v1", 0) != 0) return false;

  if (!std::getline(in, line)) return false;
  if (line != fingerprint())
  {
    std::cout << " Cache " << opt.cache << " exists but was generated with different"
              << " parameters; it will be recomputed." << std::endl;
    return false;
  }

  const int ngroups = (int) grp_type.size();

  std::string mode;
  if (!(in >> mode)) return false;

  if (mode == "rest")
  {
    rest_state.assign(ngroups, arma::vec(nvars, arma::fill::zeros));
    for (int g = 0; g < ngroups; g++)
      for (int i = 0; i < nvars; i++)
        if (!(in >> rest_state[g](i))) return false;
    traj.clear();
    snap_t.reset();
    return true;
  }
  else if (mode == "traj")
  {
    int nsnap = 0;
    if (!(in >> nsnap) || nsnap <= 0) return false;

    snap_t.set_size(nsnap);
    for (int j = 0; j < nsnap; j++)
      if (!(in >> snap_t(j))) return false;

    traj.assign(ngroups, arma::mat(nvars, nsnap, arma::fill::zeros));
    rest_state.assign(ngroups, arma::vec(nvars, arma::fill::zeros));

    for (int g = 0; g < ngroups; g++)
    {
      for (int j = 0; j < nsnap; j++)
        for (int i = 0; i < nvars; i++)
          if (!(in >> traj[g](i, j))) return false;
      rest_state[g] = traj[g].col(0);
    }
    return true;
  }

  return false;
}

void CellWarmup::save_cache() const
{
  if (opt.cache.empty() || !is_ready) return;

  std::ofstream out(opt.cache.c_str());
  if (!out.is_open())
  {
    std::cout << " *** warm up: could not write the cache "
              << opt.cache << std::endl;
    return;
  }

  out << "# cardiax cell warm-up cache v1\n";
  out << fingerprint() << "\n";
  out << std::scientific << std::setprecision(16);

  const int ngroups = (int) grp_type.size();

  if (opt.phase_from_lat && !traj.empty())
  {
    const int nsnap = (int) snap_t.n_elem;
    out << "traj\n" << nsnap << "\n";
    for (int j = 0; j < nsnap; j++) out << snap_t(j) << ((j+1)%8 ? " " : "\n");
    out << "\n";
    for (int g = 0; g < ngroups; g++)
    {
      for (int j = 0; j < nsnap; j++)
      {
        for (int i = 0; i < nvars; i++) out << traj[g](i, j) << " ";
        out << "\n";
      }
    }
  }
  else
  {
    out << "rest\n";
    for (int g = 0; g < ngroups; g++)
    {
      for (int i = 0; i < nvars; i++) out << rest_state[g](i) << " ";
      out << "\n";
    }
  }

  out.close();
  std::cout << " Limit-cycle state written to " << opt.cache << std::endl;
}
