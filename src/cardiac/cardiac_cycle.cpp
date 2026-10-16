#include <armadillo>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdexcept>
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

//! Leitor daquela TAG do pvloop do xml. Talvez tirar isso com o advento do toml?
//TODO: passar isso para o toml
struct PvloopHeader
{
  double total_time = 0.0;   //!< duracao da curva
  int    size       = 0;     //!< numero de pontos da curva
  double T_ref      = 1.0;   //!< escala de Ta dos modelos fenomenologicos
};

PvloopHeader read_pvloop_header(const std::string &mshfile)
{
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());
  if (!result)
    std::cout << "XML parsed with errors\n"
              << "Error description: " << result.description() << "\n"
              << "Error offset: " << result.offset << "\n\n";

  pugi::xml_node pv = doc.child("mesh").child("pvloop");

  PvloopHeader h;
  h.total_time = pv.attribute("total_time").as_double();
  h.size       = pv.attribute("size").as_int();
  h.T_ref      = pv.attribute("T_ref") ? pv.attribute("T_ref").as_double() : 1.0;
  return h;
}

//! -circ 1 -> Regazzoni, -pw 1 -> Piecewise, None -> legacy (Wind).
std::unique_ptr<CirculationModel> select_circulation()
{
  const int use_circ = CommandLineArgs::read("-circ", 0,
      "0D circulation in closed mesh (Regazzoni 2022)", "regazzoni");
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

CardiacElectromechanic::CardiacElectromechanic(const std::string &epmodel) : timer()
{

  // ~~~
}

CardiacElectromechanic::~CardiacElectromechanic()
{
  // ~~~
}

// ======================================================================
//  Configuracao
// ======================================================================
void CardiacElectromechanic::config(const std::string &basename)
{
  //TODO: migrar para toml algum desses arquivos. 
  std::string mshfile  = basename;
  std::string parfile  = basename;
  std::string monofile = basename;
  std::string cellmodel, odesolver;
  filename = basename;

  cout << endl << "Starting Cardiac Cycle problem" << endl;

  cellmodel = CommandLineArgs::read("-c", "Kerkoff2003");   // Kerkoff2003 como padrao
  odesolver = CommandLineArgs::read("-m", "ExplicitEuler");

  // -------------------------------------------------------------------
  //  Tempo. TRES quantidades, TRES donos:
  //
  //    periodo cardiaco   <- o modelo de circulacao (-hr/-bpm, ou a malha)
  //    duracao da rodada  <- -beats * periodo, LIMITADA por -t
  //    BCL do warm up     <- -bcl, default = periodo cardiaco
  //
  //  Antes de separar, -t acumulava as tres, e elas so eram compativeis com
  //  -beats 1: -t curto deixava o relogio da EP terminar no meio da rodada
  //  (Eikonal::advance() virava um no-op SILENCIOSO), e -t longo dava ao
  //  warm up um BCL de varios batimentos.
  // -------------------------------------------------------------------
  const PvloopHeader pv = read_pvloop_header(mshfile);
  T_ref = pv.T_ref;

  TimeBase tb;
  tb.unit      = CommandLineArgs::read("-tunit", std::string("s"),
      "Eletrophysiology time unity (s or ms)", "time"); // Segundos é padrão
  tb.solver_ms = (tb.unit == "ms") ? 1.0 : 1000.0; // 1000 é o padrão
  tb.ceiling   = CommandLineArgs::read("-t", -1.0,
      "Simulation time; with -circ 1 or -pw 1 this is optional", "time"); // tempo total de simulação
  tb.beats     = CommandLineArgs::read("-beats", 1, 
      "Number of beats", "circulation");

  // The EP step comes from -dt; with out the -dt flag, it comes from <pvloop> of the xml 
  const double dt_cli = CommandLineArgs::read("-dt", -1.0,
      "EP Timestep, in -tunit unity of measurement", "time");
      const double dt_curve = pv.total_time / (pv.size - 1);
  if (dt_cli > 0.0)
    cout << " dt from the pvloop curve (" << dt_curve
         << ") overridden by -dt " << dt_cli << endl;
  tb.dt = (dt_cli > 0.0) ? dt_cli : dt_curve;
  if (!(tb.dt > 0.0))
    throw std::runtime_error("Invalid EP time step: use -dt, or add <pvloop> tag in xml");

  circulation = select_circulation();
  circulation->configure(mshfile, tb);
  const RunPlan &plan = circulation->plan();
  cout << " Circulation model: " << circulation->name() << endl;

  const double dt = tb.dt;
  const double T  = plan.run_ms / tb.solver_ms + plan.span_margin_steps * dt;

  if (plan.period_ms > 0.0)
  {
    const double bcl_cli = CommandLineArgs::read("-bcl", -1.0,
        "BCL do warm up celular (default: periodo cardiaco)", "warmup");
    const double bcl = (bcl_cli > 0.0) ? bcl_cli : plan.period_ms / tb.solver_ms;
    ephy.set_warmup_bcl(bcl);

    cout << " Heart period = " << plan.period_ms / tb.solver_ms << " " << tb.unit
         << " | run = " << plan.beats << " beat(s)"
         << " | EP span = " << T << " " << tb.unit
         << " | warm-up BCL = " << bcl << " " << tb.unit << endl;
  }

  cout << "total time: " << T << endl;
  cout << "Mechanical time step: " << dt << endl;

  ephy.setup(monofile, cellmodel, odesolver, dt, T, 1.0, 1.0);
  ephy.update_matrix(false);
  ephy.init();

  // Coordenadas Cobiveco: opcionais. Se a malha trouxer o bloco <tm>, o tipo
  // celular (endo/mid/epi) passa a vir da coordenada transmural em vez de ser
  // uniforme. Tem de vir depois de init() -- que e onde Cells e alocado -- e
  // antes de initial_conditions(), porque Cells::init() ja usa o tipo para
  // escolher as condicoes iniciais de cada celula.
  ephy.read_cobiveco(mshfile);
  if (ephy.has_cobiveco())
    ephy.set_cell_types_from_tm(CommandLineArgs::read("-tm_endo", 0.3),
                                CommandLineArgs::read("-tm_epi",  0.7));

  // Gradiente apicobasal de IKs, ligado por default quando a malha traz <ab>.
  // -abgrad 0 desliga, para poder comparar as duas rodadas na mesma malha.
  if (CommandLineArgs::read("-abgrad", 1))
    ephy.set_apicobasal_from_ab();
  else
    cout << " Apicobasal gradient disabled by -abgrad 0" << endl;

  //TODO: conferir isso para todos os modelos (principalmente ToRORd e pro Kerkoff)
  //TODO: conferir se o código faz sentido: a ideia é ter dt diferente para mecânico e EP, e rodar várias vezes o EP para cada passo mecânico. 
  // Declara a unidade de tempo da EP. Modelos ionicos como o ToRORd sao
  // escritos em ms e precisam de passo ~0.01-0.02 ms com Euler explicito, o
  // que e impraticavel em segundos; o Kerckhoffs e escrito em segundos. Cells
  // converte o tempo para a unidade de cada modelo, entao ambos funcionam
  // sob qualquer das duas:
  //   -tunit ms  -> -t 800 -dt 0.02
  //   -tunit s   -> -t 0.8 -dt 0.001
  ephy.set_solver_time_unit_ms(tb.solver_ms);
  cout << " EP time unit: " << tb.unit << " (1 unit = " << tb.solver_ms
       << " ms); cell model " << cellmodel << " is written in "
       << (cellmodel == "Kerkoff2003" ? "s" : "ms") << endl;

  // Sanidade: T deve ser um ou poucos batimentos, ou seja, algumas centenas
  // de ms. Uma divergencia aqui significa que a unidade declarada e a
  // configuracao real discordam, e toda conversao de tempo adiante estaria
  // errada por um fator 1000.
  const double T_ms = T * tb.solver_ms;
  cout << " EP span: T = " << T << " " << tb.unit << " = " << T_ms
       << " ms, dt = " << dt << " " << tb.unit << " = " << dt * tb.solver_ms
       << " ms" << endl;
  const double T_ms_max = 20000.0 * std::max(1, plan.beats);
  if (T_ms < 50.0 || T_ms > T_ms_max)
  {
    cout << "\n *** TIME UNIT WARNING ***" << endl;
    cout << " With -tunit " << tb.unit << ", the EP span works out to " << T_ms
         << " ms, which is not a plausible number of beats." << endl;
    cout << " Either -tunit is wrong, or T came from the mesh's"
         << " <pvloop total_time=...> in a different unit." << endl;
    cout << " Pass -t and -dt explicitly in " << tb.unit << " to override.\n"
         << endl;
  }

  configure_units(cellmodel);
  configure_mef(cellmodel);

  elas.config(mshfile, parfile);
  elas.set_output_step(false);
  elas.init();

  // Dimensiona o writer uma unica vez: o numero de frames sai do plano da
  // rodada. +2: o frame 0 e escrito antes de qualquer solucao, e o ultimo
  // passo pode cair fora da divisao exata.
  out_every = CommandLineArgs::read("-prc", plan.default_out_every);
  log_every = CommandLineArgs::read("-logc", 25);
  if (out_every < 1) out_every = 1;
  if (log_every < 1) log_every = 1;
  n_frames = plan.nsteps / out_every + 2;
  cout << " Output: " << plan.nsteps << " steps, saving every " << out_every
       << " -> " << n_frames << " frame(s)" << endl;
  elas.setup_data_writer(n_frames);

  configure_fibers();

  // Init ta and vm vectors with n_points
  const int npoints = elas.get_mesh().get_n_points();
  ta.zeros(npoints);
  vm_node.zeros(npoints);

  // Dimensiona os buffers de estiramento para o estado indeformado. O
  // primeiro frame sai antes de qualquer solucao mecanica, e um campo que
  // falta no frame 0 mas aparece depois deixa a serie XDMF inconsistente.
  if (lambda_coupling)
  {
    lambda_node.ones(npoints);
    lambda_rate_node.zeros(npoints);
  }

  //TODO: não sei se deve ser mantido assim, eu acho que ou é definido por material, ou é definido globalmente... não sei se precisa dos dois. Acho redundante. 
  // Contratilidade global. Separada do ta_scale por material de proposito: o
  // primeiro e calibracao fisiologica e muda a cada rodada; o segundo e
  // atributo anatomico, fixo pela geometria.
  ta_scale_global = CommandLineArgs::read("-tascale", 1.0,
      "escala GLOBAL de contratilidade, multiplica Ta em todos os nos",
      "mecanica");

  //TODO: mesma coisa de cima... parece redundante. Mas talvez a ideia aqui seja o contrário: ao invés de ser global, ser local por nó. 
  // Mapa do ta_scale por elemento: constante no tempo, calculado uma vez.
  elas.active_scale_field(ta_scale_map);
  describe_active_tension_pipeline();

  // Read eikonal solution from file
  ephy.read_eikonal_solution(basename);

  if (CommandLineArgs::read("-dumpflags", 1,
        "print all flags used after config()", "Diagnosis"))
    CommandLineArgs::dump(cout);
}

// TODO: To achando muito estranho tudo sobre essas unidades de medida... Não está claro quais valores colocar no "padrão"
void CardiacElectromechanic::configure_units(const std::string &cellmodel)
{
  // Unidades de pressao do modelo por kPa: 1000 se a mecanica roda em Pa, 1 se
  // em kPa. Converte a tensao ativa de Land, que e dada em kPa, para a unidade
  // que o material e as condicoes de pressao usam.
  kPa_to_p = CommandLineArgs::read("-kpa2p", 1000.0);

  // A circulacao trabalha em mmHg. Essa conversao deriva da MESMA escolha de
  // unidade de pressao, em vez de ser uma flag independente: duas flags que
  // descrevem um sistema de unidades podem discordar, e um fator 1000
  // silencioso entre a tensao ativa e as pressoes das cavidades seria muito
  // dificil de ver nos resultados.
  mmHg_to_p = 0.133322387415 * kPa_to_p;
  cout << " Pressure unit: 1 kPa = " << kPa_to_p << " model units, "
       << "1 mmHg = " << mmHg_to_p << " model units" << endl;

  cell = CellProfile::for_model(cellmodel, kPa_to_p, T_ref, ephy.ms_to_solver_time());
  land.set_state_indices(cell.land_xs, cell.land_xw, cell.land_zetas, cell.land_zetaw);
}

void CardiacElectromechanic::configure_fibers()
{
  // Fiber-sheet-normal directions for mechanics
  for (int i = 0; i < elas.get_mesh().get_n_elements(); i++)
  {
    arma::vec3 f0 = ephy.get_fiber(i);
    arma::vec3 s0 = ephy.get_trans(i);
    arma::vec3 n0 = ephy.get_normal(i);
    elas.get_mesh().set_element_fiber(i, f0);
    elas.get_mesh().set_element_trans(i, s0);
    elas.get_mesh().set_element_normal(i, n0);
  }

  // Meanwhile the EP and the EPHY meshes must be the same. 
  assert(elas.get_mesh().get_n_points() == ephy.get_mesh().get_n_points());
  assert(elas.get_mesh().get_n_elements() == ephy.get_mesh().get_n_elements());
}

void CardiacElectromechanic::describe_active_tension_pipeline()
{
  if (ta_scale_map.n_elem > 0)
  {
    const double kpa2p = cell.ta_unit;
    const double smax  = ta_scale_map.max();

    cout << "\n Active tension pipeline:" << endl;
    cout << "   1. cell model              Ta [kPa]" << endl;
    cout << "   2. x " << kpa2p
         << " (kPa -> model pressure units)" << endl;
    cout << "   3. x " << ta_scale_global
         << " (-tascale, global contractility, uniform over all nodes)"
         << endl;
    cout << "        -> logged as 'Ta max ... nodal', saved as 'active_stress'"
         << endl;
    cout << "   4. x ta_scale of the element's material, in ["
         << ta_scale_map.min() << ", " << smax << "] (per element, exact,"
         << " NOT interpolated)" << endl;
    cout << "        -> logged as '... applied', saved as 'Ta_applied'" << endl;
    cout << "   net: Ta of 1 kPa in the cell model becomes up to "
         << kpa2p * ta_scale_global * smax << " model pressure units" << endl;

    if (smax > 1.5)
    {
      cout << "\n *** ATENCAO: o maior ta_scale de material e " << smax
           << "." << endl;
      cout << "     Com -tascale, os valores por material sao FRACOES de"
           << " contratilidade regional" << endl;
      cout << "     (1 = miocardio normal, 0 = passivo). Um valor > 1.5"
           << " costuma ser resto da" << endl;
      cout << "     convencao antiga, e multiplica com -tascale "
           << ta_scale_global << " dando " << ta_scale_global * smax
           << "x." << endl;
      cout << "     Se a intencao era " << smax
           << "x de contratilidade, ponha ta_scale=1 no XML e use -tascale "
           << smax << "." << endl;
    }
  }
}

void CardiacElectromechanic::solve()
{
  const RunPlan &plan = circulation->plan();
  const double ep_per_ms = ephy.ms_to_solver_time();

  cout << "\n=== " << circulation->name() << " ===" << endl;
  timer.enter("solve");

  elas.pre_solve();
  ephy.initial_conditions();
  update_active_tension();            // Ta em t = 0

  circulation->begin(*this);

  const int n_sub = cell_substeps(plan, ep_per_ms);
  ensure_ep_covers(plan, ep_per_ms);
  const double step_ep = plan.step_ms * ep_per_ms;

  if (cell.stim_amplitude != 0.0)
    cout << " LAT stimulus: amplitude = " << cell.stim_amplitude
         << ", duration = " << cell.stim_duration
         << ", period = " << plan.period_ms * ep_per_ms
         << ", offset = " << plan.stim_offset_ms * ep_per_ms
         << " (solver units)" << endl;
  else if (plan.stim_offset_ms != 0.0)
    cout << " No LAT stimulus: the cell model is driven directly by the LAT"
         << " (Kerckhoffs), so activation does NOT honour the phase offset;"
         << " consider shifting the mesh LAT by " << plan.stim_offset_ms
         << " ms." << endl;

  int frame = 0;
  write_frame(frame, plan);
  cout << " Output: saving every " << out_every << " step(s), " << n_frames
       << " frame(s) allocated" << endl;
  cout << " Log: detailed report every " << log_every << " step(s)" << endl;
  cout << "\n Legend: [ODE] cell model" << endl;

  double t_ms = 0.0;
  for (int i = 0; i < plan.nsteps; i++)
  {
    const bool log_now = (i % log_every == 0) || (i == plan.nsteps - 1);
    quiet_solve = !log_now;

    t_ms += plan.step_ms;
    const StepContext step = { i, plan.nsteps, t_ms, plan.step_ms, log_now };
    circulation->announce(step);

    update_fiber_stretch(step_ep);

    // TODO: isso daqui é preciso compreender se eu avanço um passo de tempo apenas, ou se eu avanço N tal que N = dt_mech/dt_ephy
    const bool stimulating = advance_cells(n_sub, plan, ep_per_ms);

    if (!ta.is_finite())
    {
      cout << "\n *** the active tension contains NaN/Inf at t = " << t_ms
           << " ms." << endl;
      cout << "     The cell model has diverged." <<endl; 
      break;
    }

    if (log_now)
    {
      // Log of cellular model
      cout << "  [ODE  ] " << n_sub << " sub-steps"
           << (stimulating ? " (stimulus on)" : "")
           << " -> Ta max = " << ta.max() << " nodal";
      if (ta_applied.n_elem > 0)
        cout << ", " << ta_applied.max() << " applied";
      cout << endl;
    }

    circulation->advance(step);

    if (i % out_every == 0 && frame < n_frames - 1)
      write_frame(++frame, plan);
  }

  quiet_solve = false;
  circulation->end();

  timer.leave();
  elas.timer.summary();
  timer.summary();
}

int CardiacElectromechanic::cell_substeps(const RunPlan &plan, double ep_per_ms)
{
  // O passo de acoplamento tem de ser multiplo do passo da EP; senao o
  // numero de sub-passos arredonda e a tensao ativa passa a andar numa
  // frequencia ligeiramente diferente da do modelo de circulacao.
  const double dt_ephy = ephy.get_time_parameters().get_dt();
  const double step_ep = plan.step_ms * ep_per_ms;
  const int    n_sub   = std::max(1, static_cast<int>(std::round(step_ep / dt_ephy)));
  const double mismatch = std::fabs(n_sub * dt_ephy - step_ep);

  cout << " Cell model: dt = " << dt_ephy << ", coupling step = " << step_ep
       << " (same units) -> " << n_sub << " sub-step(s) per step" << endl;
  if (mismatch > 1e-9 * step_ep)
  {
    cout << " *** WARNING: the coupling step is not a multiple of the"
         << " cell-model step." << endl;
    cout << "     " << n_sub << " x " << dt_ephy << " = " << n_sub * dt_ephy
         << " vs " << step_ep << "  (drift of " << mismatch << " per step)."
         << endl;
  }
  return n_sub;
}

void CardiacElectromechanic::ensure_ep_covers(const RunPlan &plan, double ep_per_ms)
{
  const double dt_ephy  = ephy.get_time_parameters().get_dt();
  const double total_ep = plan.run_ms * ep_per_ms;
  double stop = ephy.get_time_parameters().get_stop();

  // Se estiver curto, estica: quem manda no tempo e este laco, nao o -t. A
  // excecao abaixo so dispara se a aritmetica ficar inconsistente, e e MUITO
  // preferivel ao no-op silencioso de Eikonal::advance(), que congela as
  // celulas sem emitir nada.
  if (stop < total_ep - 0.5 * dt_ephy)
  {
    const double novo = total_ep + 2.0 * dt_ephy;
    cout << " EP span esticado de " << stop << " para " << novo
         << " para cobrir a rodada" << endl;
    ephy.time_parameters().set_stop(novo);
    stop = ephy.get_time_parameters().get_stop();
  }
  if (stop < total_ep - 0.5 * dt_ephy)
    throw std::runtime_error("EP span menor que a rodada: o modelo celular "
                             "pararia de integrar no meio, sem avisar. "
                             "Verifique -t, -hr/-pwrr, -pwinit e -beats.");
}

bool CardiacElectromechanic::advance_cells(int n_sub, const RunPlan &plan,
                                           double ep_per_ms)
{
  const double period = plan.period_ms * ep_per_ms;
  const double offset = plan.stim_offset_ms * ep_per_ms;
  const bool   stim   = (cell.stim_amplitude != 0.0);

  bool stimulating = false;
  const arma::vec &lat = ephy.get_lat();

  for (int k = 0; k < n_sub; k++)
  {
    if (stim)
    {
      // "estimulo ligado" tem de refletir se algum no foi de fato estimulado
      // neste passo, e nao apenas que o estimulo esta habilitado.
      double tk = ephy.get_time_parameters().time() - offset;
      if (period > 0.0)
        tk = std::fmod(std::fmod(tk, period) + period, period);
      if (tk >= 0.0 && arma::any(tk >= lat && tk < lat + cell.stim_duration))
        stimulating = true;

      if (offset == 0.0)
        ephy.apply_lat_stimulus(cell.stim_amplitude, cell.stim_duration, period);
      else
        ephy.apply_lat_stimulus(cell.stim_amplitude, cell.stim_duration, period, offset);
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

void CardiacElectromechanic::apply_pressures(double /*t_ms*/, double p_lv, double p_rv)
{
  if (!quiet_solve)
    cout << "Pressure: LV=" << p_lv << " RV=" << p_rv
         << " Ta: min=" << ta.min() << " max=" << ta.max() << endl;

  elas.set_pressure_Ta(kLvMarker, p_lv, kRvMarker, p_rv, ta);

  elas.set_active_stabilization(ka, lambda_raw);

  elas.solve();
  elas.reset();
}

double CardiacElectromechanic::volume_LV()    { return elas.volume_LV(); }
double CardiacElectromechanic::volume_RV()    { return elas.volume_RV(); }
double CardiacElectromechanic::volume_total() { return elas.total_volume_cavity(); }

//TODO: testar ver se está funcionando com todos os modelos celulares. Principalmente Kerkoff e Torord
void CardiacElectromechanic::update_active_tension()
{
  // Modelos fenomenologicos (Kerckhoffs) publicam a tensao ativa como valor
  // monitorado 0. Modelos ionicos (ToRORdLand) nao: o Ta deles e uma variavel
  // de estado (49), que tem de ser lida diretamente.
  if (cell.ta_var < 0)
    ephy.get_cells().get_monitored_values(0, ta);
  else
    ephy.get_cells().get_var(cell.ta_var, ta);

  // Duas escalas, com papeis distintos:
  //   cell.ta_unit    = conversao kPa -> unidade de pressao do modelo
  //   ta_scale_global = contratilidade global, -tascale
  // Ambas uniformes, aplicadas no campo NODAL. A terceira escala, a do
  // material, e por elemento e entra so na montagem.
  ta = ta * cell.ta_unit * ta_scale_global;

  // Tensao efetivamente entregue a montagem: a media nodal de Ta em cada
  // elemento vezes o ta_scale do material daquele elemento. Por elemento,
  // porque ta_scale e constante por elemento -- um no na fronteira entre
  // marcadores nao tem valor efetivo unico.
  elas.effective_active_tension(ta, ta_applied);

  update_active_stiffness();

  // So diagnostico -- os valores NAO sao modificados aqui. Uma tensao ativa
  // negativa e uma tensao compressiva na fibra, que e o que inverte
  // elementos alguns passos depois; esta linha e o aviso que precede o
  // Jacobiano negativo.
  if (lambda_rate_on && ta.n_elem > 0 && ta.min() < 0.0)
    land.report_negative_ta(ephy.get_cells(), ta, lambda_rate_node,
                            lambda_rate_max_clip);
}

void CardiacElectromechanic::update_active_stiffness()
{
  // A estabilizacao so tem serventia quando a tensao ativa de fato depende da
  // velocidade de encurtamento. Sem -lamrate, ZETAS e ZETAW ficam em repouso,
  // Ta nao tem dependencia de velocidade, e nao ha laco de realimentacao a
  // estabilizar -- entao a mecanica fica exatamente como era.
  if (!stabilize_active || !lambda_rate_on || !land.has_stiffness_states()
      || lambda_raw.n_elem != ta.n_elem)
  {
    ka.reset();
    return;
  }

  land.active_stiffness(ephy.get_cells(), lambda_raw, ka);

  // Ka sai do submodelo de Land em kPa, exatamente como Ta, e portanto
  // precisa da MESMA conversao que ta recebeu; uma escala diferente deixaria
  // um fator 5000 entre o termo de estabilizacao e a tensao que ele
  // estabiliza.
  ka = ka * cell.ta_unit;
}


// TODO: essa parte está funcionando? Já com o writer HDF5 do jeito que está? 
void CardiacElectromechanic::save_node_fields(int frame)
{
  // Both quantities live on the nodes: Cells is built with one ODE system
  // per mesh point (Eikonal::init), and the active tension is a nodal field
  // in the material (allocate_Ta(npoints)).

  // Membrane potential -- ionic models only. Kerckhoffs has no Vm, and its
  // state 0 is the contractile element length, which must not be written
  // into a field named vm.
  if (cell.vm_var >= 0)
  {
    ephy.get_cells().get_var(cell.vm_var, vm_node);
    elas.store_point_field(frame, vm_node, string("vm"));
  }

  // Active tension, in the same unit the mechanics sees: `ta` has already
  // been scaled by set_active_tension_source (kPa -> model pressure unit for
  // ToRORdLand, T_ref for Kerckhoffs). Divide by kPa_to_p when
  // post-processing if you want kPa.
  // Ta nodal: modelo celular x (kPa -> unidades) x -tascale global. Uniforme,
  // sem a escala do material. Divida por kPa_to_p e por -tascale para voltar
  // aos kPa do modelo celular.
  elas.store_point_field(frame, ta, string("active_stress"));

  // Tensao ativa EFETIVA, por elemento. "active_stress" acima e o valor
  // nodal bruto do modelo celular; "Ta_applied" ja inclui o ta_scale do
  // material e e a grandeza que aparece no residuo. Onde ta_scale = 0
  // (plugs valvares, cicatriz) este campo e identicamente zero, o que torna
  // a regiao passiva visivel no ParaView.
  if (ta_applied.n_elem > 0)
    elas.store_cell_field(frame, ta_applied, string("Ta_applied"));

  // Mapa do multiplicador. Constante no tempo, gravado a cada frame para o
  // XDMF ficar consistente e permitir sobrepor aos demais campos.
  if (ta_scale_map.n_elem > 0)
    elas.store_cell_field(frame, ta_scale_map, string("ta_scale"));

  // Fibre stretch, when the mechano-electric coupling is on. Saving it makes
  // it possible to see whether lambda is doing anything sensible before
  // trying to interpret its effect on Ta.
  if (lambda_coupling && lambda_node.n_elem == ta.n_elem)
    elas.store_point_field(frame, lambda_node, string("lambda_f"));

  // The rate is the field that actually decides whether Ta can go negative,
  // so it is worth having offline. Only written when the rate coupling is
  // on, which keeps the output identical for every existing run.
  if (lambda_coupling && lambda_rate_on
      && lambda_rate_node.n_elem == ta.n_elem)
    elas.store_point_field(frame, lambda_rate_node, string("lambda_rate"));
}

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
  // (1/ms for ToRORd), while dt_solver is in the EP solver's unit, which is
  // seconds or milliseconds depending on -tunit. Cells knows the ratio, so
  // ask it rather than reproducing the conversion here: getting this wrong
  // is a silent factor of 1000 on ZETAS/ZETAW.
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
