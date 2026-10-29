#include <armadillo>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <highfive/highfive.hpp>
#include "util/command_line_args.h"
#include "util/pugixml.hpp"
#include "cardiac_cycle.hpp"
#include "legacy_circulation.hpp"
#include "regazzoni_circulation.hpp"
#include "piecewise_circulation.hpp"

namespace
{
// Marcadores das superficies das cavidades na malha mecanica
constexpr int kLvMarker = 30;
constexpr int kRvMarker = 20;

// Unico item do <pvloop> usado fora do LegacyCirculation: a escala de Ta dos
// modelos fenomenologicos. Com o TOML vira [cell_model].ta_ref.
double read_T_ref(const std::string &mshfile)
{
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());
  if (!result)
    std::cout << "XML parsed with errors: " << result.description() << "\n";

  pugi::xml_attribute a = doc.child("mesh").child("pvloop").attribute("T_ref");
  return a ? a.as_double() : 1.0;
}

// -circ 1 -> Regazzoni, -pw 1 -> por fases, nenhum -> legado
std::unique_ptr<CirculationModel> select_circulation()
{
  const int use_circ = CommandLineArgs::read("-circ", 0,
      "0D circulation in closed loop (Regazzoni 2022)", "regazzoni");
  const int use_pw = CommandLineArgs::read("-pw", 0,
      "Piecewise pvloop (Wang et al. 2025 / Alya)", "piecewise");

  if (use_circ && use_pw)
    throw std::runtime_error("-circ 1 and -pw 1 are mutually exclusive: "
                             "both prescribe endocardial pressure.");

  if (use_circ) return std::unique_ptr<CirculationModel>(new RegazzoniCirculation());
  if (use_pw)   return std::unique_ptr<CirculationModel>(new PiecewiseCirculation());
  return std::unique_ptr<CirculationModel>(new LegacyCirculation());
}
} // namespace

CardiacElectromechanic::CardiacElectromechanic(const std::string &epmodel) : timer() {}

CardiacElectromechanic::~CardiacElectromechanic() {}

// ======================================================================
//  Configuracao
// ======================================================================
void CardiacElectromechanic::config(const std::string &basename)
{
  // O mesmo XML e malha, material e setup da EP: um caminho so (nao-const,
  // porque algumas funcoes abaixo recebem std::string&).
  std::string mesh = basename;
  std::string cellmodel = CommandLineArgs::read("-c", "Kerkoff2003");
  std::string odesolver = CommandLineArgs::read("-m", "ExplicitEuler");
  filename = basename;

  cout << endl << "Starting Cardiac Cycle problem" << endl;

  T_ref = read_T_ref(mesh);
  cell  = CellProfile::for_model(cellmodel, T_ref);
  land.set_state_indices(cell.land_xs, cell.land_xw, cell.land_zetas, cell.land_zetaw);

  // ---- tempo (ms) ----------------------------------------------------
  // A cada passo de acoplamento (-dt_mech) a EP avanca N = dt_mech/dt passos
  // (-dt). -t so LIMITA a rodada: ela e -beats periodos, e um -t curto reduz
  // os batimentos em vez de deixar o relogio da EP parar no meio.
  TimeBase tb;
  tb.ceiling_ms = CommandLineArgs::read("-t", -1.0,
      "Simulation time CEILING [ms]; with -circ 1 or -pw 1 it is optional", "time");
  tb.step_ms    = CommandLineArgs::read("-dt_mech", 1.0,
      "Coupling step: cell model -> circulation -> mechanics [ms]", "time");
  tb.beats      = CommandLineArgs::read("-beats", 1, "Number of beats", "circulation");
  tb.period_ms  = CommandLineArgs::read("-rr", -1.0,
      "Cardiac period [ms] (default: the circulation model's own, 1000)", "circulation");
  if (!(tb.step_ms > 0.0))
    throw std::runtime_error("-dt_mech must be > 0 (ms)");

  circulation = select_circulation();
  circulation->configure(mesh, tb);
  const RunPlan &plan = circulation->plan();
  cout << " Circulation model: " << circulation->name() << endl;

  // Passo da EP: -dt, ou o que o modelo de circulacao sugere (legado: um passo
  // por amostra da curva), ou o default do modelo celular
  const double dt_cli = CommandLineArgs::read("-dt", -1.0,
      "EP time step [ms]; default depends on the cell model", "time");
  const double dt = (dt_cli > 0.0)             ? dt_cli
                  : (plan.ep_dt_hint_ms > 0.0) ? plan.ep_dt_hint_ms
                                               : cell.default_dt_ms;

  // Span da EP = rodada + folga: finished() e (cur_time >= stop), entao sem
  // folga o ultimo sub-passo seria descartado
  const double T = plan.run_ms + plan.span_margin_steps * dt;

  if (plan.period_ms > 0.0)
  {
    const double bcl_cli = CommandLineArgs::read("-bcl", -1.0,
        "Cell warm-up BCL [ms] (default: heart period)", "warmup");
    const double bcl = (bcl_cli > 0.0) ? bcl_cli : plan.period_ms;
    ephy.set_warmup_bcl(bcl);
    cout << " Heart period = " << plan.period_ms << " ms | run = " << plan.beats
         << " beat(s) | EP span = " << T << " ms | warm-up BCL = " << bcl << " ms" << endl;
  }
  cout << " EP step = " << dt << " ms | coupling step = " << plan.step_ms
       << " ms | " << plan.nsteps << " coupling step(s)" << endl;

  ephy.setup(mesh, cellmodel, odesolver, dt, T, 1.0, 1.0);
  ephy.update_matrix(false);
  ephy.init();

  // Cobiveco (opcional): com o bloco <tm> na malha, o tipo celular vem da
  // coordenada transmural. Tem de vir depois de init() (onde Cells e alocado)
  // e antes de initial_conditions() (Cells::init() ja usa o tipo).
  ephy.read_cobiveco(mesh);
  if (ephy.has_cobiveco())
    ephy.set_cell_types_from_tm(CommandLineArgs::read("-tm_endo", 0.3),
                                CommandLineArgs::read("-tm_epi",  0.7));

  // Gradiente apicobasal de IKs: ligado por default quando a malha traz <ab>
  if (CommandLineArgs::read("-abgrad", 1))
    ephy.set_apicobasal_from_ab();
  else
    cout << " Apicobasal gradient disabled by -abgrad 0" << endl;

  // A EP roda em ms (como Monodomain::init()); Cells converte para a unidade
  // nativa de cada modelo. O default de Cells e 1000 (s): sem esta linha o
  // ToRORd rodaria 1000x lento, sem aviso.
  ephy.set_solver_time_unit_ms(1.0);

  configure_mef(cellmodel);

  elas.config(mesh, mesh);
  elas.set_output_step(false);
  elas.save_checkpoint(false);   // o checkpoint e do cardiac_cycle, nao por incremento
  elas.init();

  // Frames: +2 porque o frame 0 sai antes de qualquer solucao
  out_every = std::max(1, CommandLineArgs::read("-prc", plan.default_out_every));
  log_every = std::max(1, CommandLineArgs::read("-logc", 25));
  n_frames  = plan.nsteps / out_every + 2;
  cout << " Output: " << plan.nsteps << " steps, saving every " << out_every
       << " -> " << n_frames << " frame(s)" << endl;
  elas.setup_data_writer(n_frames);

  configure_fibers();

  // Ta e um campo NODAL: Cells tem um sistema de EDOs por ponto da malha
  const int npoints = elas.get_mesh().get_n_points();
  ta.zeros(npoints);
  vm_node.zeros(npoints);

  // Buffers de estiramento no estado indeformado: o frame 0 sai antes de
  // qualquer solucao mecanica, e um campo que falta nele deixa o XDMF inconsistente
  if (lambda_coupling)
  {
    lambda_node.ones(npoints);
    lambda_rate_node.zeros(npoints);
  }

  const double tascale = CommandLineArgs::read("-tascale", 1.0,
      "Global contractility scale (multiplies Ta at every node)", "mechanics");
  ta_gain = cell.ta_unit * tascale;

  // Ta prescrita: a tabela ja esta em Pa, entao o ganho e so o -tascale
  const std::string ta_curve_file = CommandLineArgs::read("-ta_curve", "",
      "Prescribed Ta(t) table (t [ms], Ta [Pa]), repeated every period; replaces"
      " the cell model's Ta", "mechanics");
  if (!ta_curve_file.empty())
  {
    ta_curve.load(ta_curve_file, plan.period_ms);
    ta_gain = tascale;
  }

  elas.active_scale_field(ta_scale_map);     // ta_scale do material, por elemento
  describe_active_tension_pipeline(tascale);

  ephy.read_eikonal_solution(basename);      // LAT por no, em ms

  checkpoint_interval = CommandLineArgs::read("-save_state", -1.0,
      "Save a checkpoint every N ms (checkpoint_t_<t>ms.h5 in the CWD)", "checkpoint");
  restore_file = CommandLineArgs::read("-restore", std::string(""),
      "Restart from a checkpoint written by -save_state (same flags and mesh)",
      "checkpoint");

  if (CommandLineArgs::read("-dumpflags", 1,
        "print all flags used after config()", "Diagnosis"))
    CommandLineArgs::dump(cout);
}

void CardiacElectromechanic::configure_fibers()
{
  // Direcoes fibra-folha-normal para a mecanica
  for (int i = 0; i < elas.get_mesh().get_n_elements(); i++)
  {
    arma::vec3 f0 = ephy.get_fiber(i);
    arma::vec3 s0 = ephy.get_trans(i);
    arma::vec3 n0 = ephy.get_normal(i);
    elas.get_mesh().set_element_fiber(i, f0);
    elas.get_mesh().set_element_trans(i, s0);
    elas.get_mesh().set_element_normal(i, n0);
  }

  // As malhas mecanica e eletrica tem de ser a mesma
  assert(elas.get_mesh().get_n_points() == ephy.get_mesh().get_n_points());
  assert(elas.get_mesh().get_n_elements() == ephy.get_mesh().get_n_elements());
}

void CardiacElectromechanic::describe_active_tension_pipeline(double tascale)
{
  if (ta_scale_map.n_elem == 0) return;

  const double smax = ta_scale_map.max();
  cout << "\n Active tension pipeline:\n";
  if (ta_curve.active())
    cout << "   Ta (-ta_curve, Pa) x " << tascale << " (-tascale) at every node  -> 'active_stress'\n";
  else
    cout << "   Ta (cell model) x " << cell.ta_unit << " (-> Pa) x " << tascale
         << " (-tascale) = x" << ta_gain << " at every node  -> 'active_stress'\n";
  cout
       << "   x ta_scale of the element's material, in [" << ta_scale_map.min()
       << ", " << smax << "] (per element)  -> 'Ta_applied'" << endl;

  // O ta_scale do material e uma FRACAO regional (1 = miocardio, 0 = passivo);
  // a calibracao global e o -tascale
  if (smax > 1.5)
    cout << " *** WARNING: largest material ta_scale is " << smax << ". If you meant "
         << smax << "x contractility, set ta_scale=1 in the XML and use -tascale "
         << smax << "." << endl;
}

// ======================================================================
//  O laco: generico para qualquer CirculationModel
// ======================================================================

void CardiacElectromechanic::solve()
{
  const RunPlan &plan = circulation->plan();

  cout << "\n=== " << circulation->name() << " ===" << endl;
  timer.enter("solve");

  // pre_solve() monta as matrizes; initial_conditions() reinicia as celulas.
  // Pular qualquer um deixa o sistema linear sem montagem (PETSc reason -11).
  elas.pre_solve();
  ephy.initial_conditions();
  update_active_tension();            // Ta em t = 0

  circulation->begin(*this);

  const int n_sub = cell_substeps(plan);

  if (cell.stim_amplitude != 0.0)
    cout << " LAT stimulus: amplitude = " << cell.stim_amplitude << ", duration = "
         << cell.stim_duration << " ms, period = " << plan.period_ms
         << " ms, offset = " << plan.stim_offset_ms << " ms" << endl;
  else if (plan.stim_offset_ms != 0.0)
    cout << " No LAT stimulus: the cell model is driven by the LAT itself"
         << " (Kerckhoffs), so activation does NOT honour the phase offset;"
         << " consider shifting the mesh LAT by " << plan.stim_offset_ms << " ms." << endl;

  // Com -restore, o estado inteiro (celulas, mecanica, circulacao, laco) vem
  // do checkpoint e sobrescreve o que foi inicializado acima
  int frame = 0;
  int first = 0;
  double t_ms = 0.0;
  if (!restore_file.empty())
    first = restore_checkpoint(restore_file, t_ms, frame);
  else
    write_frame(frame, plan);

  // Checkpoint a cada N passos de acoplamento
  const int checkpoint_rate = (checkpoint_interval > 0.0)
      ? std::max(1, static_cast<int>(std::round(checkpoint_interval / plan.step_ms))) : -1;

  cout << " Output every " << out_every << " step(s), log every " << log_every
       << " step(s)\n\n Legend: [ODE] cell model" << endl;

  for (int i = first; i < plan.nsteps; i++)
  {
    const bool log_now = (i % log_every == 0) || (i == plan.nsteps - 1);
    quiet_solve = !log_now;

    t_ms += plan.step_ms;
    const StepContext step = { i, plan.nsteps, t_ms, plan.step_ms, log_now };
    circulation->announce(step);

    // 0. Realimentacao mecano-eletrica. lambda_f vem da ULTIMA solucao mecanica
    //    convergida: uma vez por passo, antes das celulas, e nunca dentro de uma
    //    inversao de pressao (o estiramento de uma pressao de teste rejeitada
    //    nao tem significado).
    update_fiber_stretch(plan.step_ms);

    // 1. Celulas: N = dt_mech/dt passos da EP. Ta fica congelado enquanto a
    //    circulacao itera.
    const bool stimulating = advance_cells(n_sub, plan);

    // Um modelo ionico rigido com passo grande demais diverge primeiro em Vm;
    // Ta segue decaindo suave por um tempo, e a rodada parece viva com a EP morta
    if (!ta.is_finite())
    {
      cout << "\n *** the active tension contains NaN/Inf at t = " << t_ms
           << " ms: the cell model has diverged. Reduce -dt." << endl;
      break;
    }

    if (log_now)
    {
      // 'ta' e do modelo celular x ganho escalar; 'applied' inclui tambem o
      // ta_scale do material, por elemento. Nao sao proporcionais.
      elas.effective_active_tension(ta, ta_applied);
      cout << "  [ODE  ] " << n_sub << " sub-steps" << (stimulating ? " (stimulus on)" : "")
           << " -> Ta max = " << ta.max() << " nodal";
      if (ta_applied.n_elem > 0) cout << ", " << ta_applied.max() << " applied";
      cout << endl;
    }

    if (i == 0 && !ta_curve.active() && arma::max(arma::abs(ta)) == 0.0)
    {
      cout << " *** WARNING: the active tension is zero at the start of the run.\n"
           << "     With Kerckhoffs, Ta is non-zero only for lat <= t <= lat + 483 ms,"
           << " where lat comes from the per-node LAT in the mesh or <pvloop passive_time=...>."
           << endl;
      if (plan.period_ms > 0.0)
        cout << "     Check that passive_time is smaller than the heart period ("
             << plan.period_ms << " ms), or the ventricle will fill and never eject." << endl;
    }

    // 2. Circulacao + mecanica
    circulation->advance(step);

    // 3. Saida
    if (i % out_every == 0 && frame < n_frames - 1)
      write_frame(++frame, plan);

    // 4. Checkpoint
    if (checkpoint_rate > 0 && (i + 1) % checkpoint_rate == 0)
      save_checkpoint(i, t_ms, frame);
  }

  quiet_solve = false;
  circulation->end();

  timer.leave();
  elas.timer.summary();
  timer.summary();
}

int CardiacElectromechanic::cell_substeps(const RunPlan &plan)
{
  // N = dt_mech/dt, ARREDONDADO (o electromech trunca: 0.3/0.1 = 2.9999999999999996
  // perderia um passo da EP). Se dt_mech nao for multiplo de dt, a tensao ativa
  // andaria numa frequencia ligeiramente diferente da da circulacao.
  const double dt_ephy  = ephy.get_time_parameters().get_dt();
  const int    n_sub    = std::max(1, static_cast<int>(std::round(plan.step_ms / dt_ephy)));
  const double mismatch = std::fabs(n_sub * dt_ephy - plan.step_ms);

  cout << " Cell model: dt = " << dt_ephy << " ms, coupling step = " << plan.step_ms
       << " ms -> " << n_sub << " sub-step(s) per step" << endl;
  if (mismatch > 1e-9 * plan.step_ms)
    cout << " *** WARNING: -dt_mech is not a multiple of -dt: " << n_sub << " x " << dt_ephy
         << " = " << n_sub * dt_ephy << " vs " << plan.step_ms << " (drift of "
         << mismatch << " ms per step)." << endl;

  // O relogio da EP tem de cobrir a rodada, senao Eikonal::advance() vira um
  // no-op silencioso. config() ja o dimensiona pelo plano; isto so dispara se
  // alguem mudar uma das pontas.
  if (ephy.get_time_parameters().get_stop() < plan.run_ms - 0.5 * dt_ephy)
    throw std::runtime_error("EP span shorter than the run: the cell model would "
                             "stop integrating mid-run, silently.");
  return n_sub;
}

bool CardiacElectromechanic::advance_cells(int n_sub, const RunPlan &plan)
{
  const bool stim = (cell.stim_amplitude != 0.0);
  bool stimulating = false;
  const arma::vec &lat = ephy.get_lat();

  for (int k = 0; k < n_sub; k++)
  {
    if (stim)
    {
      // "estimulo ligado" = algum no foi de fato estimulado neste passo
      double tk = ephy.get_time_parameters().time() - plan.stim_offset_ms;
      if (plan.period_ms > 0.0)
        tk = std::fmod(std::fmod(tk, plan.period_ms) + plan.period_ms, plan.period_ms);
      if (tk >= 0.0 && arma::any(tk >= lat && tk < lat + cell.stim_duration))
        stimulating = true;

      if (plan.stim_offset_ms == 0.0)
        ephy.apply_lat_stimulus(cell.stim_amplitude, cell.stim_duration, plan.period_ms);
      else
        ephy.apply_lat_stimulus(cell.stim_amplitude, cell.stim_duration,
                                plan.period_ms, plan.stim_offset_ms);
    }
    ephy.advance();
  }

  update_active_tension();
  return stimulating;
}

void CardiacElectromechanic::write_frame(int frame, const RunPlan &plan)
{
  elas.output_vtk(0, frame);
  if (plan.store_stress) elas.storeStress(frame);
  save_node_fields(frame);
}

// ======================================================================
//  Checkpoint (-save_state / -restore)
// ======================================================================

void CardiacElectromechanic::save_checkpoint(int step, double t_ms, int frame)
{
  timer.enter("Save Checkpoint");
  const TimeParameters &tip = ephy.get_time_parameters();
  const Cells &cells = ephy.get_cells();

  // 1. EP + mecanica (x, fext0), como no electromech
  arma::vec vm(elas.get_mesh().get_n_points());
  cells.get_var(0, vm);
  elas.save_coupled_checkpoint(tip.it(), tip.time(), vm.memptr(),
                               cells.get_state_vars(), cells.get_ode_size());

  // 2. O resto, no mesmo arquivo. O nome e o de WriterHDF5::write_coupled_checkpoint.
  char time_buf[64];
  std::snprintf(time_buf, sizeof(time_buf), "%.2f", tip.time());
  const std::string file = "checkpoint_t_" + std::string(time_buf) + "ms.h5";
  cout << " Saving checkpoint at t = " << t_ms << " ms (step " << step + 1 << "/"
       << circulation->plan().nsteps << ") -> " << file << endl;

  HighFive::File h5(file, HighFive::File::ReadWrite);

  // Mecanica: o que save_coupled_checkpoint nao grava. eL entra na pressao da
  // proxima solucao; F e lido por update_fiber_stretch() ANTES dela.
  const std::vector<arma::mat33*> &vecF = elas.get_vec_F();
  std::vector<double> F(9 * vecF.size());
  for (size_t k = 0; k < vecF.size(); k++)
    std::copy(vecF[k]->memptr(), vecF[k]->memptr() + 9, F.begin() + 9 * k);
  HighFive::Group gm = h5.getGroup("/mechanics");
  gm.createDataSet("lagrange_multipliers",
                   arma::conv_to<std::vector<double>>::from(elas.lagrange_multipliers()));
  gm.createDataSet("F", F);

  // Laco e realimentacao mecano-eletrica
  HighFive::Group gc = h5.createGroup("/cardiac_cycle");
  gc.createAttribute("next_step", step + 1);
  gc.createAttribute("t_ms", t_ms);
  gc.createAttribute("frame", frame);
  gc.createAttribute("circulation", std::string(circulation->name()));
  gc.createAttribute("lambda_step_count", lambda_step_count);
  gc.createAttribute("lambda_prev_valid", lambda_prev_valid ? 1 : 0);
  const std::pair<const char*, const arma::vec*> mef[] = {
    { "lambda_node", &lambda_node }, { "lambda_raw", &lambda_raw },
    { "lambda_prev", &lambda_prev }, { "lambda_rate_node", &lambda_rate_node } };
  for (const auto &m : mef)
    if (m.second->n_elem > 0)
      gc.createDataSet(m.first, arma::conv_to<std::vector<double>>::from(*m.second));

  // Circulacao
  CirculationModel::State st;
  circulation->save_state(st);
  HighFive::Group gq = h5.createGroup("/circulation");
  for (const auto &kv : st)
    gq.createDataSet(kv.first, kv.second);

  timer.leave();
}

int CardiacElectromechanic::restore_checkpoint(const std::string &file, double &t_ms, int &frame)
{
  cout << "\n Restoring checkpoint " << file << endl;

  HighFive::File h5(file, HighFive::File::ReadOnly);
  if (!h5.exist("/cardiac_cycle"))
    throw std::runtime_error("Checkpoint " + file + " was not written by cardiac_cycle"
                             " (no /cardiac_cycle group)");
  HighFive::Group gc = h5.getGroup("/cardiac_cycle");

  std::string model;
  gc.getAttribute("circulation").read(model);
  if (model != circulation->name())
    throw std::runtime_error("Checkpoint written with circulation model '" + model
                             + "', but this run uses '" + circulation->name() + "'");

  const Cells &cells = ephy.get_cells();
  const std::vector<size_t> sv = h5.getDataSet("/ep/state_variables").getSpace().getDimensions();
  if (sv.size() != 2 || sv[0] * sv[1] != cells.get_size() || sv[1] != cells.get_ode_size())
    throw std::runtime_error("Checkpoint cell states do not match this mesh/cell model");

  // 1. EP + mecanica, como no electromech
  int ep_step;
  double ep_time;
  arma::vec vm(sv[0]);
  elas.read_ep_checkpoint(file, vm.memptr(), cells.get_state_vars(), ep_step, ep_time);
  elas.restore_checkpoint(file);
  elas.lc.reset();
  const_cast<TimeParameters &>(ephy.get_time_parameters()).restore_state(ep_step, ep_time);

  HighFive::Group gm = h5.getGroup("/mechanics");
  std::vector<double> eL, F;
  gm.getDataSet("lagrange_multipliers").read(eL);
  gm.getDataSet("F").read(F);
  const std::vector<arma::mat33*> &vecF = elas.get_vec_F();
  if (eL.size() != elas.lagrange_multipliers().n_elem || F.size() != 9 * vecF.size())
    throw std::runtime_error("Checkpoint mechanics does not match this mesh");
  elas.lagrange_multipliers() = arma::vec(eL);
  for (size_t k = 0; k < vecF.size(); k++)
    std::copy(F.begin() + 9 * k, F.begin() + 9 * (k + 1), vecF[k]->memptr());

  // 2. Laco e realimentacao mecano-eletrica
  int next_step, prev_valid;
  gc.getAttribute("next_step").read(next_step);
  gc.getAttribute("t_ms").read(t_ms);
  gc.getAttribute("frame").read(frame);
  gc.getAttribute("lambda_step_count").read(lambda_step_count);
  gc.getAttribute("lambda_prev_valid").read(prev_valid);
  lambda_prev_valid = (prev_valid != 0);
  const std::pair<const char*, arma::vec*> mef[] = {
    { "lambda_node", &lambda_node }, { "lambda_raw", &lambda_raw },
    { "lambda_prev", &lambda_prev }, { "lambda_rate_node", &lambda_rate_node } };
  for (const auto &m : mef)
    if (gc.exist(m.first))
    {
      std::vector<double> v;
      gc.getDataSet(m.first).read(v);
      *m.second = arma::vec(v);
    }

  // 3. Circulacao (depois de begin(), que ela sobrescreve)
  CirculationModel::State st;
  HighFive::Group gq = h5.getGroup("/circulation");
  for (const std::string &name : gq.listObjectNames())
    gq.getDataSet(name).read(st[name]);
  circulation->restore_state(st);

  // Ta e Ka do estado restaurado (as mesmas do fim do passo salvo)
  update_active_tension();

  cout << " Restarting at t = " << t_ms << " ms, step " << next_step + 1 << "/"
       << circulation->plan().nsteps << " (EP t = " << ep_time << " ms), frame "
       << frame << endl;
  return next_step;
}

// ======================================================================
//  A porta para os modelos de circulacao (CavityMechanics)
// ======================================================================

void CardiacElectromechanic::apply_pressures(double /*t_ms*/, double p_lv, double p_rv)
{
  if (!quiet_solve)
    cout << "Pressure: LV=" << p_lv << " RV=" << p_rv << " Pa  Ta: min=" << ta.min()
         << " max=" << ta.max() << endl;

  elas.set_pressure_Ta(kLvMarker, p_lv, kRvMarker, p_rv, ta);

  // Estabilizacao da tensao ativa dependente da velocidade. lambda_raw, nao
  // lambda_node: o incremento (lambda_atual - lambda_prev) tem de diferenciar
  // dois valores da mesma regra, e o material calcula lambda_atual da
  // deformacao bruta. Contra o campo CORTADO haveria um salto quando um no
  // entra ou sai de [-lammin, -lammax], multiplicado por uma rigidez ~10^6.
  // Vazios, os dois vetores restauram o comportamento sem estabilizacao.
  elas.set_active_stabilization(ka, lambda_raw);

  elas.solve();
  elas.reset();
}

// Os volumes saem de UpdatedLagrangian em mL (malha em um, pressao em Pa).
double CardiacElectromechanic::volume_LV()    { return elas.volume_LV(); }
double CardiacElectromechanic::volume_RV()    { return elas.volume_RV(); }
double CardiacElectromechanic::volume_total() { return elas.total_volume_cavity(); }

// ======================================================================
//  Tensao ativa
// ======================================================================

void CardiacElectromechanic::update_active_tension()
{
  // Quem sabe de onde vem Ta e o CellProfile: fenomenologicos (Kerckhoffs) o
  // publicam como valor monitorado 0; ionicos (ToRORdLand), no estado 49.
  // Com -ta_curve, a Ta e a da tabela, igual em todos os nos.
  if (ta_curve.active())
    ta.fill(ta_curve(ephy.get_time_parameters().time()));
  else if (cell.ta_var < 0)
    ephy.get_cells().get_monitored_values(0, ta);
  else
    ephy.get_cells().get_var(cell.ta_var, ta);

  // Um ganho escalar (conversao x -tascale), no campo NODAL. A escala do
  // material, por elemento, entra so na montagem.
  ta = ta * ta_gain;

  update_active_stiffness();

  // So diagnostico. Ta negativa e tensao compressiva na fibra, o que inverte
  // elementos alguns passos depois; este e o aviso que precede o Jacobiano negativo.
  if (lambda_rate_on && ta.n_elem > 0 && ta.min() < 0.0)
    land.report_negative_ta(ephy.get_cells(), ta, lambda_rate_node, lambda_rate_max_clip);
}

void CardiacElectromechanic::update_active_stiffness()
{
  // So ha o que estabilizar quando Ta depende da velocidade de encurtamento,
  // ou seja, com -lamrate. Sem ele a mecanica fica exatamente como era.
  if (!stabilize_active || !lambda_rate_on || !land.has_stiffness_states()
      || lambda_raw.n_elem != ta.n_elem)
  {
    ka.reset();
    return;
  }

  land.active_stiffness(ephy.get_cells(), lambda_raw, ka);

  // Ka e dTa/dlambda e a mecanica recebe Ta com o ganho INTEIRO: Ka precisa do
  // mesmo, senao o termo de estabilizacao sai de proporcao com a tensao que
  // estabiliza.
  ka = ka * ta_gain;
}

void CardiacElectromechanic::save_node_fields(int frame)
{
  // Vm so para modelos ionicos: o estado 0 do Kerckhoffs e o comprimento do
  // elemento contratil e nao pode ir para um campo "vm"
  if (cell.vm_var >= 0)
  {
    ephy.get_cells().get_var(cell.vm_var, vm_node);
    elas.store_point_field(frame, vm_node, string("vm"));
  }

  // Ta nodal: modelo celular x ganho, sem a escala do material
  elas.store_point_field(frame, ta, string("active_stress"));

  // Ta EFETIVA por elemento (inclui o ta_scale do material, e e a que entra no
  // residuo). Zero onde ta_scale = 0 (plugs, cicatriz): torna a regiao passiva
  // visivel no ParaView.
  elas.effective_active_tension(ta, ta_applied);
  if (ta_applied.n_elem > 0)
    elas.store_cell_field(frame, ta_applied, string("Ta_applied"));

  // Mapa do ta_scale: constante, gravado a cada frame para o XDMF ficar consistente
  if (ta_scale_map.n_elem > 0)
    elas.store_cell_field(frame, ta_scale_map, string("ta_scale"));

  if (lambda_coupling && lambda_node.n_elem == ta.n_elem)
    elas.store_point_field(frame, lambda_node, string("lambda_f"));

  // So com o acoplamento de taxa ligado, para a saida das rodadas existentes
  // continuar identica
  if (lambda_coupling && lambda_rate_on && lambda_rate_node.n_elem == ta.n_elem)
    elas.store_point_field(frame, lambda_rate_node, string("lambda_rate"));
}

// ======================================================================
//  Realimentacao mecano-eletrica (estiramento de fibra). Sem alteracao.
// ======================================================================

void CardiacElectromechanic::configure_mef(const std::string &cellmodel)
{
  // ---- mechano-electric feedback: fibre stretch ----------------------
  // Off by default, so every existing run reproduces exactly as before.
  //   -lam 1      feeds lambda_f = sqrt(I4f) to the cell model
  //   -lamrate 1  additionally feeds d(lambda_f)/dt
  //   -lamclip a b  overrides the guard rails on the value handed over
  //   -lamratemax r caps |d(lambda)/dt| at r, in the CELL MODEL's time unit
  //                 (1/ms for ToRORd). 0 disables the cap.
  lambda_coupling = (CommandLineArgs::read("-lam", 0) != 0);
  lambda_rate_on  = (CommandLineArgs::read("-lamrate", 0) != 0);
  lambda_min_clip = CommandLineArgs::read("-lammin", 0.7);
  lambda_max_clip = CommandLineArgs::read("-lammax", 1.3);
  lambda_rate_max_clip = CommandLineArgs::read("-lamratemax", 0.003);
  lambda_rate_delay    = CommandLineArgs::read("-lamratedelay", 0);
  stabilize_active     = (CommandLineArgs::read("-lamstab", 1) != 0);

  if (lambda_coupling)
  {
    // Only the Land submodel reads the stretch. Silently accepting -lam on a
    // model that ignores it would look like the coupling was active.
    if (cellmodel != "ToRORdLand")
    {
      cout << " *** -lam 1 requested, but cell model '" << cellmodel
           << "' does not use stretch; coupling disabled." << endl;
      lambda_coupling = false;
    }
    else
    {
      cout << " Mechano-electric feedback ON: lambda_f = sqrt(I4f)"
           << " per node -> ToRORd-Land" << endl;
      cout << "   lambda clipped to [" << lambda_min_clip << ", "
           << lambda_max_clip << "]" << endl;
      cout << "   d(lambda)/dt: "
           << (lambda_rate_on ? "on" : "off (isometric contraction as far as"
                                       " the cell is concerned)")
           << endl;
      if (lambda_rate_on)
      {
        // Ta = Lfac*(Tref/dr)*[ (ZETAS+1)*XS + ZETAW*XW ], and in the
        // relaxed limit ZETAS -> 249.1*rate, ZETAW -> 24.6*rate. Ta changes
        // sign when
        //     rate < -XS / (249.1*XS + 24.6*XW).
        // With XS >> XW that reduces to rate < -1/249.1 = -0.004, the
        // ZETAS-driven case the old message described. But in diastole and
        // early systole XS is SMALL while XW stays around 0.15, and then the
        // ZETAW*XW term dominates: at XS = 0.01, XW = 0.15 the threshold is
        // already -0.0016. The single constant cap cannot express that, so
        // print both regimes instead of only the favourable one.
        if (lambda_rate_max_clip > 0.0)
        {
          cout << "   |d(lambda)/dt| limitado a " << lambda_rate_max_clip
               << " (in the cell model's time unit)" << endl;
          cout << "   Ta changes sign for rate < -XS/(249.1*XS + 24.6*XW):"
               << endl;
          cout << "     XS=0.20, XW=0.15 (systole)  -> threshold -0.0037"
               << (lambda_rate_max_clip < 0.0037 ? "  [covered]"
                                                 : "  [NOT covered]") << endl;
          cout << "     XS=0.01, XW=0.15 (diastole) -> threshold -0.0016"
               << (lambda_rate_max_clip < 0.0016 ? "  [covered]"
                                                 : "  [NOT covered]") << endl;
        }
        else
          cout << "   *** |d(lambda)/dt| UNBOUNDED (-lamratemax 0): the active"
               << " tension may turn negative and invert elements" << endl;

        if (lambda_rate_delay > 0)
          cout << "   d(lambda)/dt held at zero for the first "
               << lambda_rate_delay << " step(s) (passive inflation)" << endl;

        cout << "   Regazzoni-Quarteroni stabilisation: "
             << (stabilize_active ? "ON (Ta -> Ta + Ka*dlambda)"
                                  : "*** OFF (-lamstab 0): plain segregated"
                                    " scheme, unstable whenever Ka > Kp")
             << endl;
      }
    }
  }
  else
    cout << " Mechano-electric feedback off (use -lam 1 to enable)"
         << endl;
}

void CardiacElectromechanic::update_fiber_stretch(double dt_solver)
{
  if (!lambda_coupling) return;

  // The Land equations expect d(lambda)/dt in the model's OWN time unit
  // (1/ms for ToRORd), while dt_solver is in the EP solver's unit (ms).
  // Cells knows the ratio to the model's unit (1 for ToRORd, 1e-3 for
  // Kerckhoffs, written in seconds), so ask it rather than reproducing the
  // conversion here: getting this wrong is a silent factor of 1000 on
  // ZETAS/ZETAW.
  const double dt_cell = dt_solver * ephy.get_cells().time_factor();

  // lambda_f = sqrt(I4f) per element, averaged onto the nodes, since the
  // cell models are nodal.
  elas.fiber_stretch_nodes(lambda_node);

  // NaN/Inf are not a matter of taste: they would poison both the value fed
  // to the cell model and the backward difference below. Replace them by the
  // undeformed state first, so that everything downstream sees finite data.
  int n_nonfinite = 0;
  for (arma::uword i = 0; i < lambda_node.n_elem; i++)
    if (!std::isfinite(lambda_node(i))) { lambda_node(i) = 1.0; n_nonfinite++; }

  // Keep the UNCLIPPED field: the rate must be a difference of two values
  // produced by the same rule. Differencing a clipped field manufactures a
  // step whenever a node enters or leaves the clip, and that step is then
  // divided by dt -- a large spurious lambda_rate for purely bookkeeping
  // reasons.
  lambda_raw = lambda_node;

  // Guard rails. A diverging Newton or an inverted element can produce a
  // lambda far outside the physiological range, and the Land submodel would
  // then be evaluated where it was never fitted. Clipping keeps a bad
  // mechanical step from destroying the cell state as well.
  int n_clipped = n_nonfinite;
  for (arma::uword i = 0; i < lambda_node.n_elem; i++)
  {
    double & l = lambda_node(i);
    if (l < lambda_min_clip)      { l = lambda_min_clip; n_clipped++; }
    else if (l > lambda_max_clip) { l = lambda_max_clip; n_clipped++; }
  }

  if (n_clipped > 0 && !quiet_solve)
    cout << "  [MEF  ] " << n_clipped << "/" << lambda_node.n_elem
         << " node(s) had lambda clipped to ["
         << lambda_min_clip << ", " << lambda_max_clip << "]" << endl;

  // d(lambda)/dt by backward difference, in the cell model's time unit.
  //
  // Held at zero for the first lambda_rate_delay calls: the passive inflation
  // moves lambda from 1 to ~1.1 over the first few steps, and differencing
  // that gives a rate an order of magnitude past the sign-change threshold
  // while nothing physiological is happening yet. lambda_prev keeps being
  // updated during the delay, so the first real rate is a proper one-step
  // difference and not a jump accumulated over the whole warm-up.
  lambda_step_count++;
  const bool rate_active = lambda_rate_on
                        && (lambda_step_count > lambda_rate_delay);

  if (lambda_rate_on && !rate_active && !quiet_solve)
    cout << "  [MEF  ] d(lambda)/dt suppressed (step "
         << lambda_step_count << "/" << lambda_rate_delay << ")" << endl;

  int n_rate_clipped = 0;
  if (rate_active && dt_cell > 0.0 && lambda_prev_valid
      && lambda_prev.n_elem == lambda_raw.n_elem)
  {
    lambda_rate_node = (lambda_raw - lambda_prev) / dt_cell;

    // Cap on |d(lambda)/dt|. Two reasons, one physical and one numerical.
    //
    // Physical: the Land distortion equation gives ZETAS_ss = (A/cds)*rate
    // = 249*rate with the stock parameters, so rate = -0.004 (1/ms) puts
    // ZETAS at -1, where Ta = Lfac*(Tref/dr)*((ZETAS+1)*XS + ZETAW*XW)
    // changes sign. A negative Ta is a compressive fibre stress and inverts
    // elements. That rate corresponds to 4 lengths/s, i.e. the model's own
    // Vmax -- so the cap enforces a bound the model already implies rather
    // than introducing a new calibrated constant.
    //
    // Numerical: a backward difference amplifies noise as 1/dt. Solver
    // tolerances, the element->node averaging and a stiff pressure inversion
    // all leave O(1e-3) jitter in lambda, which at dt = 0.1 ms is already a
    // rate of 1e-2 -- 2.5x past the sign change. This is why refining the
    // time step does not help: it makes the difference quotient worse, not
    // better.
    if (lambda_rate_max_clip > 0.0)
    {
      const double r_max = lambda_rate_max_clip;
      for (arma::uword i = 0; i < lambda_rate_node.n_elem; i++)
      {
        double & r = lambda_rate_node(i);
        if (r > r_max)       { r =  r_max; n_rate_clipped++; }
        else if (r < -r_max) { r = -r_max; n_rate_clipped++; }
      }
    }
  }
  else
    lambda_rate_node.zeros(lambda_raw.n_elem);

  // The history is the unclipped field, for the reason given above.
  lambda_prev = lambda_raw;
  lambda_prev_valid = true;

  ephy.set_fiber_stretch(lambda_node, lambda_rate_node);

  if (!quiet_solve)
  {
    cout << "  [MEF  ] lambda_f: min = " << lambda_node.min()
         << ", mean = " << arma::mean(lambda_node)
         << ", max = " << lambda_node.max() << endl;

    if (rate_active)
    {
      // Both distortion states, not just ZETAS. They relax towards
      //   ZETAS_ss = 249.1*rate  (tau_s = 24.9 ms, slow to build up)
      //   ZETAW_ss =  24.6*rate  (tau_w =  2.46 ms, essentially immediate)
      // so ZETAS is the larger number but ZETAW is the one already at its
      // steady value within a single dtc. Printing both makes it visible
      // which term is loaded before Ta actually changes sign.
      const double r_abs = arma::max(arma::abs(lambda_rate_node));
      cout << "  [MEF  ] d(lambda)/dt: min = " << lambda_rate_node.min()
           << ", max = " << lambda_rate_node.max() << endl;
      cout << "  [MEF  ]   worst case: ZETAS_ss = " << -249.128 * r_abs
           << ", ZETAW_ss = " << -24.639 * r_abs << endl;
      if (n_rate_clipped > 0)
        cout << "  [MEF  ] " << n_rate_clipped << "/"
             << lambda_rate_node.n_elem << " node(s) with |d(lambda)/dt|"
             << " capped at " << lambda_rate_max_clip << endl;
    }
  }
}