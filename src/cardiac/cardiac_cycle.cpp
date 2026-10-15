#include <armadillo>
#include <cmath>
#include <iomanip>
#include <ios>
#include <stdexcept>
#include <algorithm>
#include "util/command_line_args.h"
#include "odes/torord_land.hpp"   // land_active_stiffness()
#include "cardiac_cycle.hpp"
#include "util/pugixml.hpp"


CardiacElectromechanic::CardiacElectromechanic(const std::string &epmodel) : dt_mech(1.0), timer(), Ta0(0.0), dTa0(0.0)
{
  // ~~
}

CardiacElectromechanic::~CardiacElectromechanic()
{
  // ~~
}

void CardiacElectromechanic::config(const string &basename)
{
  double T, dt;
  string parfile = basename;
  string mshfile = basename;
  string monofile = basename;
  string cellmodel, odesolver;
  filename = basename;

  cout << endl
       << "Starting coupled electromechanics problem" << endl;

  T = CommandLineArgs::read("-t", 1.0);
  dt = CommandLineArgs::read("-dt", 0.1);
  cellmodel = CommandLineArgs::read("-c", "Kerkoff2003"); // using kerkoff2003 as a standard value
  odesolver = CommandLineArgs::read("-m", "ExplicitEuler");

  double pr = 1.0;
  double pa = 1.0;

  ifstream inp;

  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());

  if (!result)
  {
    std::cout << "XML parsed with errors\n";
    std::cout << "Error description: " << result.description() << "\n";
    std::cout << "Error offset: " << result.offset;
    std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
  }

  pugi::xml_node fnodes = doc.child("mesh").child("pvloop");
  T = fnodes.attribute("total_time").as_double();
  int size = fnodes.attribute("size").as_int();
  c_art = fnodes.attribute("C_art").as_double();
  c_ven = fnodes.attribute("C_ven").as_double();
  part = fnodes.attribute("p_art").as_double();
  pven = fnodes.attribute("p_ven").as_double();
  Rmv = fnodes.attribute("R_mv").as_double();
  Rven = fnodes.attribute("R_ven").as_double();
  Rper = fnodes.attribute("R_per").as_double();
  Rao = fnodes.attribute("R_ao").as_double();
  int biv = fnodes.attribute("biv").as_int();

  V_art_zero = fnodes.attribute("V_art_zero").as_double();
  V_ven_zero = fnodes.attribute("V_ven_zero").as_double();

  E_es_LA = fnodes.attribute("E_es_LA").as_double();
  A_LA = fnodes.attribute("A_LA").as_double();
  B_LA = fnodes.attribute("B_LA").as_double();
  Tmax = fnodes.attribute("T_max").as_double();
  tau = fnodes.attribute("tau").as_double();

  P_o = fnodes.attribute("P_o").as_double();
  stroke_volume = fnodes.attribute("stroke_volume").as_double();
  T_ref = (fnodes.attribute("T_ref")) ? fnodes.attribute("T_ref").as_double() : 1.0;

  double press = 0, press2 = 0;
  for (pugi::xml_node node = fnodes.child("step"); node; node = node.next_sibling("step"))
  {
    double t = node.attribute("time").as_double();
    if (biv)
    {
      press = node.attribute("pressure_lv").as_double();
      press2 = node.attribute("pressure_rv").as_double();
    }
    else
    {
      press = node.attribute("pressure").as_double();
    }

    double ta = node.attribute("active_tension").as_double() * T_ref;
    curr_time.push_back(t * 1000.0);

    p_lv.push_back(press);
    p_rv.push_back(press2);

    cout << "Pressure: " << press << " " << press2 << " Ta: " << ta << endl;
  }

  
  // The pvloop curve in the mesh file dictates T and dt. That is right for
  // the legacy curve-driven path, but the closed-loop circulation path does
  // not use the curve at all and needs its own time discretisation --
  // notably a much smaller dt for stiff ionic cell models. So an explicitly
  // given -t / -dt wins over the values implied by the curve.
  double dt_from_curve = T / (size - 1);
  double T_cli  = CommandLineArgs::read("-t",  -1.0,
      "TETO de tempo de simulacao; com -circ 1 ou -pw 1 e opcional", "tempo");
  double dt_cli = CommandLineArgs::read("-dt", -1.0,
      "passo de tempo da EP, na unidade de -tunit", "tempo");

  // dt primeiro: a margem do span da EP e expressa em passos.
  if (dt_cli > 0.0)
  {
    cout << " dt from the pvloop curve (" << dt_from_curve
         << ") overridden by -dt " << dt_cli << endl;
    dt = dt_cli;
  }
  else
    dt = dt_from_curve;

  // -------------------------------------------------------------------
  // TRES quantidades, TRES donos
  //
  //   periodo cardiaco   <- -hr / -bpm
  //   duracao da rodada  <- -beats * periodo, LIMITADA por -t
  //   BCL do warm up     <- -bcl, default = periodo cardiaco
  //
  // Antes de separar, -t acumulava as tres e elas so eram mutuamente
  // compativeis com -beats 1:
  //   -t curto  -> o relogio da EP terminava no meio da rodada e
  //                Eikonal::advance() virava um no-op SILENCIOSO;
  //   -t longo  -> o warm up recebia um BCL de varios batimentos e as
  //                condicoes iniciais do ciclo limite eram de outra
  //                frequencia.
  // -------------------------------------------------------------------
  const int use_circ_early = CommandLineArgs::read("-circ", 0,
      "circulacao 0D em malha fechada (Regazzoni 2022)", "circulacao");

  // PV loop POR FASES (Wang et al. 2025). Terceiro caminho, ao lado do
  // legado (curva tabelada no <pvloop>) e do -circ. Mutuamente exclusivo
  // com -circ: os dois prescrevem a condicao de contorno endocardica, e
  // aceitar ambos silenciosamente deixaria um deles sem efeito visivel.
  const int use_pw_early = CommandLineArgs::read("-pw", 0,
      "PV loop por fases (Wang et al. 2025 / Alya)", "pvloop por fases");

  if (use_circ_early && use_pw_early)
    throw std::runtime_error("-circ 1 e -pw 1 sao mutuamente exclusivos: "
                             "ambos prescrevem a pressao endocardica. "
                             "Escolha um dos dois caminhos.");

  // Leituras puras, sem efeito colateral: set_solver_time_unit_ms() continua
  // acontecendo mais abaixo, no lugar de sempre.
  const string tunit_early     = CommandLineArgs::read("-tunit", string("s"),
      "unidade de tempo da EP (s ou ms)", "tempo");
  const double solver_ms_early = (tunit_early == "ms") ? 1.0 : 1000.0;
  const double units_per_sec   = 1000.0 / solver_ms_early;

  if (use_circ_early)
  {
    double hr  = CommandLineArgs::read("-hr",  -1.0,
        "frequencia cardiaca em Hz (default 1.0 = 60 bpm)", "circulacao");
    double bpm = CommandLineArgs::read("-bpm", -1.0,
        "frequencia cardiaca em bpm (alternativa a -hr)", "circulacao");
    if (bpm > 0.0) hr = bpm / 60.0;
    if (hr  > 0.0) circ.set_heart_rate(hr);

    const double tact = CommandLineArgs::read("-tact", -1.0,
        "instante da ativacao eletrica ventricular, em s", "circulacao");
    if (tact >= 0.0) circ.set_activation_time(tact);

    const double emd = CommandLineArgs::read("-emd", -1.0,
        "atraso eletromecanico estimulo -> inicio de Ta, em s", "circulacao");
    if (emd >= 0.0) circ.set_emd(emd);

    cout << "\n" << circ.describe_timing() << endl;

    // Periodo cardiaco na unidade de tempo do solver.
    const double RR = circ.period() * units_per_sec;

    circ_num_beats = CommandLineArgs::read("-beats", 1,
        "numero de batimentos", "circulacao");
    if (circ_num_beats < 1) circ_num_beats = 1;

    // -t e um TETO, nunca um alvo. Se nao cobrir os batimentos pedidos,
    // reduz o numero de batimentos -- em vez de deixar o relogio da EP
    // terminar no meio de um deles, que era o modo de falha silencioso.
    if (T_cli > 0.0)
    {
      const int fits = static_cast<int>(std::floor(T_cli / RR + 1e-9));
      if (fits < circ_num_beats)
      {
        cout << " *** -t " << T_cli << " " << tunit_early << " cobre apenas "
             << fits << " batimento(s) de " << RR << " " << tunit_early
             << "; -beats " << circ_num_beats << " reduzido para "
             << std::max(1, fits) << endl;
        circ_num_beats = std::max(1, fits);
      }
    }

    // O span da EP cobre a rodada inteira mais dois passos de margem: o
    // teste e finished() == (cur_time >= stop), entao sem margem o ultimo
    // sub-passo seria descartado.
    T = circ_num_beats * RR + 2.0 * dt;

    // BCL do warm up: UM batimento, nunca a rodada inteira.
    const double bcl_cli = CommandLineArgs::read("-bcl", -1.0,
        "BCL do warm up celular (default: periodo cardiaco)", "warmup");
    const double bcl     = (bcl_cli > 0.0) ? bcl_cli : RR;
    ephy.set_warmup_bcl(bcl);

    cout << " Heart period = " << RR << " " << tunit_early
         << " | run = " << circ_num_beats << " beat(s)"
         << " | EP span = " << T << " " << tunit_early
         << " | warm-up BCL = " << bcl << " " << tunit_early << endl;
  }
  else if (use_pw_early)
  {
    // ---------------------------------------------------------------
    //  PV loop por fases: o controlador trabalha em ms, o solver na
    //  unidade de -tunit. TODA conversao acontece aqui, uma vez.
    // ---------------------------------------------------------------
    use_piecewise = true;

    pw.configure(mshfile);
    cout << "\n" << pw.describe() << endl;

    pw_dt_ms = CommandLineArgs::read("-dtp", 1.0,
        "passo do controlador do PV loop por fases, em ms",
        "pvloop por fases");
    if (pw_dt_ms <= 0.0) pw_dt_ms = 1.0;

    pw_num_beats = CommandLineArgs::read("-beats", 1,
        "numero de batimentos", "circulacao");
    if (pw_num_beats < 1) pw_num_beats = 1;

    // ms -> unidade de tempo do solver. Com -tunit ms o fator e 1; com
    // -tunit s e 1/1000.
    const double ms_to_units = 1.0 / solver_ms_early;

    // -t e um TETO, como no caminho -circ. Se nao cobrir a rodada pedida,
    // reduz o numero de batimentos em vez de deixar o relogio da EP parar
    // no meio -- Eikonal::advance() vira um no-op silencioso nesse caso.
    if (T_cli > 0.0)
    {
      const double t_init_units = pw.timing().t_init * ms_to_units;
      const double RR_units     = pw.timing().RR     * ms_to_units;
      const int fits = static_cast<int>(
          std::floor((T_cli - t_init_units) / RR_units + 1e-9));
      if (fits < pw_num_beats)
      {
        cout << " *** -t " << T_cli << " " << tunit_early << " cobre apenas "
             << fits << " batimento(s) apos a fase 0 de "
             << t_init_units << " " << tunit_early << "; -beats "
             << pw_num_beats << " reduzido para " << std::max(1, fits)
             << endl;
        pw_num_beats = std::max(1, fits);
      }
    }

    // Span da EP: rodada inteira mais dois passos de margem, pelo mesmo
    // motivo do caminho -circ (finished() testa cur_time >= stop).
    T = pw.total_time(pw_num_beats) * ms_to_units + 2.0 * dt;

    // BCL do warm up: UM batimento, nunca a rodada inteira.
    const double bcl_cli = CommandLineArgs::read("-bcl", -1.0,
        "BCL do warm up celular (default: periodo cardiaco)", "warmup");
    const double bcl = (bcl_cli > 0.0) ? bcl_cli
                                       : pw.timing().RR * ms_to_units;
    ephy.set_warmup_bcl(bcl);

    // O passo do controlador tem de ser um multiplo do passo da EP; senao
    // o numero de sub-passos arredonda e a tensao ativa passa a andar numa
    // frequencia ligeiramente diferente da do controlador.
    const double dt_pw_units = pw_dt_ms * ms_to_units;
    const double ratio = dt_pw_units / dt;
    const int    n_sub = std::max(1, static_cast<int>(std::round(ratio)));
    if (std::fabs(n_sub * dt - dt_pw_units) > 1e-9 * dt_pw_units)
    {
      cout << " *** ATENCAO: -dtp " << pw_dt_ms << " ms = " << dt_pw_units
           << " " << tunit_early << " nao e multiplo de -dt " << dt
           << "; " << n_sub << " x " << dt << " = " << n_sub * dt << endl;
    }

    cout << " Heart period = " << pw.timing().RR << " ms"
         << " | fase 0 = " << pw.timing().t_init << " ms"
         << " | run = " << pw_num_beats << " batimento(s)"
         << " | EP span = " << T << " " << tunit_early
         << " | warm-up BCL = " << bcl << " " << tunit_early << endl;
    cout << " Controlador: dt = " << pw_dt_ms << " ms -> " << n_sub
         << " sub-passo(s) de EP por passo" << endl;
  }
  else if (T_cli > 0.0 && T_cli != T)
  {
    // Caminho legado: -t continua sendo a duracao da simulacao, e o BCL do
    // warm up continua saindo dela (warmup_bcl fica em 0).
    cout << " total_time from the mesh (" << T << ") overridden by -t "
         << T_cli << endl;
    T = T_cli;
  }

  dt_mech = dt;
  cout << "size: " << size << endl;
  cout << "total time: " << T << endl;
  cout << "Mechanical time step: " << dt_mech << endl;

  ephy.setup(monofile, cellmodel, odesolver, dt, T, pr, pa);
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

  // Declare the time unit the EP problem runs in. Ionic models such as
  // ToRORd are written in ms and need a step of ~0.01-0.02 ms with explicit
  // Euler, which is impractical to express in seconds; the phenomenological
  // Kerckhoffs model is written in seconds. Cells converts time into each
  // model's own unit, so both work under either setting.
  //   -tunit ms  -> pass -t 800 -dt 0.02
  //   -tunit s   -> pass -t 0.8 -dt 0.001
  string tunit = CommandLineArgs::read("-tunit", "s");
  double solver_ms = (tunit == "ms") ? 1.0 : 1000.0;
  ephy.set_solver_time_unit_ms(solver_ms);
  // EP time units per SECOND -- the reciprocal of solver_ms, which is ms per
  // EP time unit. With tunit=ms: 1000 EP units per second. With tunit=s: 1.
  ephy_time_per_second = 1000.0 / solver_ms;
  cout << " EP time unit: " << tunit << " (1 unit = " << solver_ms
       << " ms); cell model " << cellmodel << " is written in "
       << (cellmodel == "Kerkoff2003" ? "s" : "ms") << endl;

  // Sanity check: T is one or a few beats, so in the declared unit it should
  // be a few hundred (ms) or a fraction of one (s). A mismatch here means the
  // declared unit and the actual configuration disagree, and every time
  // conversion downstream would be off by a factor of 1000.
  double T_ms = T * solver_ms;
  cout << " EP span: T = " << T << " " << tunit << " = " << T_ms
       << " ms, dt = " << dt << " " << tunit << " = " << dt * solver_ms
       << " ms" << endl;
  // O teto escala com o numero de batimentos: com o span cobrindo a rodada
  // inteira, 20 batimentos de 1 s sao 20000 ms legitimos.
  const double T_ms_max =
      20000.0 * std::max(1, std::max(circ_num_beats, pw_num_beats));
  if (T_ms < 50.0 || T_ms > T_ms_max)
  {
    cout << "\n *** TIME UNIT WARNING ***" << endl;
    cout << " With -tunit " << tunit << ", the EP span works out to " << T_ms
         << " ms, which is not a plausible number of beats." << endl;
    cout << " Either -tunit is wrong, or T came from the mesh's"
         << " <pvloop total_time=...> in a different unit." << endl;
    cout << " Pass -t and -dt explicitly in " << tunit << " to override.\n"
         << endl;
  }

  // How to obtain the active tension, and whether the cells need a stimulus.
  //   Kerckhoffs: phenomenological, Ta is monitored value 0 and the model is
  //     triggered by the activation time itself -> no stimulus current.
  //   ToRORdLand: ionic, Ta is state variable 49 and the cell must actually
  //     depolarise -> a stimulus is applied at each node once its LAT is
  //     reached. The x5 factor on Ta follows the fenicsx-pulse demo.
  // Model pressure units per kPa: 1000 if the mechanics runs in Pa, 1 if in
  // kPa. It converts the Land active tension, which is expressed in kPa,
  // into whatever unit the material and the pressure BCs use.
  kPa_to_p = CommandLineArgs::read("-kpa2p", 1000.0);

  // The circulation works in mmHg. Derive that conversion from the same
  // choice of pressure unit instead of taking it as an independent flag:
  // two flags describing one unit system can disagree, and a silent factor
  // of 1000 between the active tension and the cavity pressures would be
  // very hard to spot in the results.
  mmHg_to_p = 0.133322387415 * kPa_to_p;
  cout << " Pressure unit: 1 kPa = " << kPa_to_p << " model units, "
       << "1 mmHg = " << mmHg_to_p << " model units" << endl;

  if (cellmodel == "ToRORdLand")
  {
    // Land gives Ta in kPa. Any further gain now comes from the per-material
    // ta_scale attribute in the mesh XML instead of a constant hidden here, so
    // that a single place controls it. NOTE: ta_scale = 5.0 in the XML
    // reproduces the factor that used to sit on this line.
    set_active_tension_source(49, kPa_to_p);
    set_vm_source(0);                     // Vm is state variable 0
    // Land distortion states, for the negative-Ta diagnostic only. Indices
    // match the reads at the top of torord_land.cpp:
    //   43 = XS, 44 = XW, 47 = ZETAS, 48 = ZETAW
    set_land_state_indices(43, 44, 47, 48);
    double amp = CommandLineArgs::read("-stimamp", -53.0);
    double dur = CommandLineArgs::read("-stimdur", 2.0 * ephy.ms_to_solver_time());
    set_lat_stimulus(amp, dur);
    cout << " Active tension: state variable 49 (kPa) x " << kPa_to_p
         << " model pressure units per kPa"
         << " (further gain from the per-material ta_scale)" << endl;
    cout << " LAT stimulus: amplitude " << amp << ", duration " << dur
         << " (solver time units)" << endl;
  }
  else
  {
    set_active_tension_source(-1, 0.0);   // monitored 0, scaled by T_ref
    set_vm_source(-1);                    // no membrane potential to save
    set_lat_stimulus(0.0, 0.002);         // no stimulus needed
  }

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

  elas.config(mshfile, parfile);
  elas.set_output_step(false);
  elas.init();

  // Size the data writer here, once. For the closed-loop circulation path
  // the number of frames follows from the run parameters; the length of the
  // pressure curve (curr_time) is irrelevant there because that curve is
  // not used at all.
  // Mesma decisao tomada no topo de config(); nao reler a linha de comando.
  const bool use_circ = (use_circ_early != 0);
  circ_out_every = CommandLineArgs::read("-prc", 8);
  circ_log_every = CommandLineArgs::read("-logc", 25);
  if (circ_log_every < 1) circ_log_every = 1;
  if (circ_out_every < 1) circ_out_every = 1;

  if (use_circ)
  {
    // circ_num_beats ja foi resolvido contra o teto -t la em cima. Reler
    // -beats aqui reintroduziria a divergencia justamente no caso em que o
    // teto cortou.
    double dtc = CommandLineArgs::read("-dtc", 1.0e-3,
        "passo de tempo da circulacao, em s", "circulacao");
    int nsteps = static_cast<int>(std::round(circ.period() / dtc));
    circ_n_frames = circ_num_beats * nsteps / circ_out_every + 1;
    cout << " Circulation output: " << circ_num_beats * nsteps
         << " steps, saving every " << circ_out_every << " -> "
         << circ_n_frames << " frame(s)" << endl;
    elas.setup_data_writer(circ_n_frames);
  }
  else if (use_piecewise)
  {
    // pw_num_beats ja foi resolvido contra o teto -t na secao de tempos.
    // Reler -beats aqui reintroduziria a divergencia justamente no caso em
    // que o teto cortou.
    pw_out_every = circ_out_every;
    pw_log_every = circ_log_every;

    const int nsteps_total = static_cast<int>(
        std::round(pw.total_time(pw_num_beats) / pw_dt_ms));
    // +2: o frame 0 e escrito antes de qualquer solucao, e o ultimo passo
    // pode cair fora da divisao exata.
    pw_n_frames = nsteps_total / pw_out_every + 2;

    cout << " Piecewise output: " << nsteps_total << " steps, saving every "
         << pw_out_every << " -> " << pw_n_frames << " frame(s)" << endl;
    elas.setup_data_writer(pw_n_frames);
  }
  else
    elas.setup_data_writer(curr_time.size());

  // Configure fibers for mechanical problem
  for (int i = 0; i < elas.get_mesh().get_n_elements(); i++)
  {
    arma::vec3 f0 = ephy.get_fiber(i);
    arma::vec3 s0 = ephy.get_trans(i);
    arma::vec3 n0 = ephy.get_normal(i);
    elas.get_mesh().set_element_fiber(i, f0);
    elas.get_mesh().set_element_trans(i, s0);
    elas.get_mesh().set_element_normal(i, n0);
  }

  // Configure fiber-sheet-normal directions
  for (int i = 0; i < ephy.get_mesh().get_n_elements(); i++)
  {
    arma::mat33 *R = new arma::mat33();
    arma::mat33 *R0 = new arma::mat33();

    R->col(0) = ephy.get_fiber(i);
    R->col(1) = ephy.get_trans(i);
    R->col(2) = ephy.get_normal(i);

    R0->col(0) = ephy.get_fiber(i);
    R0->col(1) = ephy.get_trans(i);
    R0->col(2) = ephy.get_normal(i);

    vec_fib.push_back(R);
    vec_fib0.push_back(R0);
  }
  int neln = ephy.get_mesh().get_nen();
  int nelem = elas.get_mesh().get_n_elements();
  for (int i = 0; i < nelem; i++)
    for (int j = 0; j < neln; j++)
      vec_stress.push_back(new arma::mat33());

  // Check if mechanical and electrical meshes are the same
  assert(elas.get_mesh().get_n_points() == ephy.get_mesh().get_n_points());
  assert(elas.get_mesh().get_n_elements() == ephy.get_mesh().get_n_elements());

  // Ta is a NODAL field: Cells holds one ODE system per mesh point
  // (Eikonal::init), get_var/get_monitored_values fill n_points entries, and
  // the material indexes it by point number (allocate_Ta(npoints),
  // Ta(pnums[i]) in updated_lagrangian). Sizing it by n_elements happened to
  // work on tetrahedral meshes only because there are more elements than
  // nodes; on hexahedra n_points > n_elements and it wrote past the end.
  int npoints = elas.get_mesh().get_n_points();
  ta.zeros(npoints);
  vm_node.zeros(npoints);

  // Pre-size the stretch buffers to the undeformed state. The first output
  // frame is written before any mechanical solve, and a field that is absent
  // in frame 0 but present afterwards makes the XDMF series inconsistent.
  if (lambda_coupling)
  {
    lambda_node.ones(npoints);
    lambda_rate_node.zeros(npoints);
  }

  // Contratilidade global. Separada do ta_scale por material de proposito: o
  // primeiro e calibracao fisiologica e muda a cada rodada; o segundo e
  // atributo anatomico, fixo pela geometria.
  ta_scale_global = CommandLineArgs::read("-tascale", 1.0,
      "escala GLOBAL de contratilidade, multiplica Ta em todos os nos",
      "mecanica");

  // Mapa do ta_scale por elemento: constante no tempo, calculado uma vez.
  elas.active_scale_field(ta_scale_map);

  if (ta_scale_map.n_elem > 0)
  {
    const double kpa2p = (ta_scale != 0.0) ? ta_scale : T_ref;
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

    // Guarda contra escala dupla. Com o -tascale global, os valores do XML
    // passam a ser FRACOES de contratilidade regional: 1 = miocardio normal,
    // 0 = plug valvar, 0.2 = cicatriz. Um valor grande ali quase sempre e
    // resto da convencao antiga, em que o XML carregava a magnitude.
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

  //"Solve" the eikonal, to discover the lat in each element
  ephy.read_eikonal_solution(basename);

  // Registro de procedencia da rodada: toda flag lida ate aqui, com o valor
  // resolvido e a marca [cli] ou [default]. Ligado por default de proposito
  // -- o modo de falha caro nao e desconhecer uma flag, e nao perceber que
  // ela ficou no default (-t ausente, -amgx caindo no fallback, etc).
  if (CommandLineArgs::read("-dumpflags", 1,
        "imprime o registro de flags no fim de config()", "diagnostico"))
    CommandLineArgs::dump(cout);
}


void CardiacElectromechanic::Solve_System(double tt, double pressure, double pressure2)
{
  int size = elas.get_mesh().get_n_elements();
  int index = static_cast<int>(tt*1000); 

  P0 = pressure;
  
  if (!quiet_solve)
    cout << "Pressure: LV=" << pressure << " RV=" << pressure2
         << " Ta: min=" << ta.min() << " max=" << ta.max() << endl;


  elas.set_pressure_Ta(30, pressure, 20, pressure2, ta);

  // Stabilisation of the velocity-dependent active tension. lambda_node is
  // lambda at the last converged mechanics solve, i.e. exactly the lambda^(k)
  // the increment is measured from: update_fiber_stretch() runs at the top of
  // the step, before the cells and before this solve. Both vectors are empty
  // when the stabilisation does not apply, which restores the old behaviour.
  // lambda_raw, not lambda_node: the increment (lambda_current - lambda_prev)
  // must difference two values produced by the same rule, and the material
  // computes lambda_current from the raw deformation. Differencing against
  // the CLIPPED field would manufacture a jump the moment a node enters or
  // leaves [-lammin, -lammax], and Ka multiplies that jump by a stiffness of
  // order 10^6.
  elas.set_active_stabilization(ka, lambda_raw);

  elas.solve();
  elas.reset();

}


void CardiacElectromechanic::solve()
{
  cout << "Solving coupled electromechanical problem" << endl;

  int nstep;
  int i = 0, ii = 0, itempo = 0, store = 1;
  int size = elas.get_mesh().get_n_points();
  int nelem = elas.get_mesh().get_n_elements();
  double dt = CommandLineArgs::read("-dt", 0.1);
  bool phase1 = true;
  bool phase3 = false;
  double VED = 0.0;

  ofstream pv_file;
  string pvfilename = filename.c_str() + string("_pvloop.txt");
  pv_file.open(pvfilename.c_str());

  arma::vec vm(size), u_field(3 * size);
  vm.zeros();

  cout << "\nTime step of mechanics: " << dt_mech << endl;
  TimeParameters tip(ephy.get_time_parameters());
  elas.pre_solve();
  ephy.initial_conditions(); 

  p_lv.push_back(0.0);
  p_rv.push_back(0.0);

  cout << "Solve 0 " << endl;
  Solve_System(0, p_lv[0], p_rv[0]);
  volume.push_back(elas.volume_LV());
  volume_rv.push_back(elas.volume_RV());

  timePoints.push_back(0.0);
  activeStressCurve.push_back(Ta0);

  cout << "Volume: " << volume[i] << endl;
  cout << "Pressure: " << p_lv[i] << " " << p_rv[i] << endl;

  double DT = dt_mech * 1000;
  double Vf_0, Vf_1, V_LP, r_est, V_LP0 = volume[0], V_LA, V_LA0 = 10.0, V_LA_zero = 10.0, V_art, V_ven, V_art0, V_ven0, p_0, p_1, pi = 3.14159265358979323846;
  double qmv, qao, qven, qper;

  p_ven.push_back(pven);
  p_art.push_back(part);
  p_LA.push_back(0.0);

  V_art0 = p_art[0] * c_art + V_art_zero;
  V_ven0 = p_ven[0] * c_ven + V_ven_zero;

  V_LP = V_LP0;
  V_art = V_art0;
  V_ven = V_ven0;
  V_LA = V_LA0;

  double total_time = 0;

  for (int k = 0; k < n_cycles; k++)
  {
    total_time += tip.time();
    tip.reset();

    itempo = 0;  
    while (!tip.finished())
    {
      tip.increase_time();
      cout << "\nTime: " << tip.time() << endl;

      i += 1;
      itempo = itempo + 1;

      double err = 1.;
      p_0 = p_lv[i - 1];

      if (i >= 5)
      {
        double dpn = p_lv[i - 1] - p_lv[i - 2];
        double dpn1 = p_lv[i - 2] - p_lv[i - 3];
        double dpn2 = p_lv[i - 3] - p_lv[i - 4];
        double dpn3 = p_lv[i - 4] - p_lv[i - 5];
        p_1 = p_0 + (DT / 24.) * (55. * dpn / DT - 59. * dpn1 / DT + 37. * dpn2 / DT - 9. * dpn3 / DT);
      }
      else
      {
        float DP = p_ven[0];
        float DV = 251. - 199.;
        float Q = DV / 200.;
        p_1 = p_0 + (DP / DV) * Q * DT;
      }

      // Mechano-electric feedback: lambda_f from the last converged solve,
      // read before the cells advance.
      update_fiber_stretch(tip.get_dt());

      //Update the active stress value
      ephy.advance();
      ephy.get_cells().get_monitored_values(0, ta);
      ta = ta * T_ref; 
                  
      Solve_System(tip.time(), p_0, p_rv[i - 1]); 
      Vf_0 = elas.total_volume_cavity();

      int iterations = 0;

      double p_ven1 = p_ven[i - 1];
      double p_art1 = p_art[i - 1];
      double p_LA1 = p_LA[i - 1];

      while (err > 0.001)
      {
        iterations += 1;
        cout << "\n Now doing Newton iteration number " << iterations << " for pressure update. \n"
             << endl;

        cout << "Solve 2 " << endl;
        Solve_System(tip.time(), p_1, p_rv[i]);
        Vf_1 = elas.total_volume_cavity();
        double C = (p_1 - p_0) / (Vf_1 - Vf_0);

        qmv = (p_LA1 >= p_1) ? (0.5 / Rmv) * (p_LA1 + p_LA[i - 1] - p_1 - p_lv[i - 1]) : 0.0;
        qao = (p_1 >= p_art1) ? (0.5 / Rao) * (p_1 + p_lv[i - 1] - p_art1 - p_art[i - 1]) : 0.0;
        qper = (0.5 / Rper) * (p_art1 + p_art[i - 1] - p_ven1 - p_ven[i - 1]);
        qven = (0.5 / Rven) * (p_ven1 + p_ven[i - 1] - p_LA1 - p_LA[i - 1]);

        V_LA = V_LA0 + (qven - qmv) * DT;
        V_LP = V_LP0 + (qmv - qao) * DT;
        V_art = V_art0 + (qao - qper) * DT;
        V_ven = V_ven0 + (qper - qven) * DT;

        p_art1 = (V_art - V_art_zero) / c_art;
        p_ven1 = (V_ven - V_ven_zero) / c_ven;

        double pES_LA = E_es_LA * (V_LA - V_LA_zero);
        double pED_LA = A_LA * (exp(B_LA * (V_LA - V_LA_zero)) - 1.0);
        double e_func = (curr_time.at(itempo) >= 0 && curr_time.at(itempo) <= (3. / 2.) * Tmax) ? 0.5 * (sin((pi / Tmax) * curr_time.at(itempo) - pi / 2.) + 1.0) : 0.5 * exp(-(curr_time.at(itempo) - (3. / 2.) * Tmax) / tau);

        p_LA1 = e_func * pES_LA + (1.0 - e_func) * pED_LA;

        r_est = Vf_1 - V_LP;

        p_0 = p_1;
        Vf_0 = Vf_1;
        p_1 = p_1 - C * r_est;

        err = fabs((Vf_1 - V_LP) / V_LP);
        cout << "\n Error after " << iterations << " iterations = " << err << "\n"
             << endl;
      }

      p_art.push_back(p_art1);
      p_ven.push_back(p_ven1);
      V_LP0 = V_LP;
      V_art0 = V_art;
      V_ven0 = V_ven;
      p_LA.push_back(p_LA1);
      V_LA0 = V_LA;

      cout << "\n  PRESSURE UPDATE CONVERGED IN " << iterations << " iterations \n\n"
           << endl;
      cout << "  Pressure is now " << p_1 << "\n"
           << endl;
      p_lv[i] = p_1;

      volume.push_back(elas.volume_LV());
      volume_rv.push_back(elas.volume_RV());

      timePoints.push_back(tip.time() + total_time);
      activeStressCurve.push_back(Ta0);

      cout << "Volume: " << volume[i] << endl;
      cout << "Pressure: " << p_lv[i] << endl;

      pv_file << tip.time() + total_time << " " << p_lv[i] << " " << volume[i] << " "
              << Ta0 << " "
              << p_art[i] << " " << p_ven[i] << " " << p_LA[i] << " "
              << V_art0 << " " << V_ven0 << " " << V_LA0 << " " << qao << " " << qmv << " " << qven << " " << qper
              << " " << p_rv[i] << " " << volume_rv[i] << endl;

      cout << itempo << " " << tip.time() + total_time << " " << curr_time.at(itempo) << " " << p_lv[i] << " " << volume[i] << " "
           << Ta0 << " " << p_art[i] << " " << p_ven[i] << " " << p_LA[i] << " " << V_art0 << " " << V_ven0 << " " << V_LA0 << endl;

      if (k == 0)
      {
        ii++;
        cout << "XDMF saving... " << "Step: " << ii << endl;
        elas.output_vtk(0, ii);
        elas.storeStress(ii);
        save_node_fields(ii);
      }
    }
    
  }

  pv_file.close();

  // string stress_filename = filename + "_active_stress.txt";
  // saveActiveStressToFile(stress_filename);

  elas.timer.summary();
  timer.summary();
}


// ======================================================================
//  Closed-loop 3D-0D coupling  (Regazzoni et al. 2022)
// ======================================================================
//
//  The 0D model owns the cavity volumes and asks the 3D model, at every
//  circulation step, "which pressures reproduce these volumes?".
//
//  fenicsx-pulse answers that directly because it constrains the volume with
//  a Lagrange multiplier, and the multiplier IS the cavity pressure. Cardiax
//  applies pressure as a Neumann condition and computes the volume, so the
//  map runs the other way and has to be inverted numerically.
//
//  The inversion is a 2x2 Newton on
//        F(p_lv, p_rv) = ( V_LV(p) - V_LV_target ,
//                          V_RV(p) - V_RV_target )
//  The off-diagonal terms matter: the septum is shared, so raising the LV
//  pressure changes the RV volume. Treating the chambers independently
//  converges poorly for exactly that reason.
//
//  The Jacobian is reused between steps and only rebuilt when convergence
//  degrades, which keeps the cost near one 3D solve per Newton iteration.

void CardiacElectromechanic::update_active_tension()
{
  // Phenomenological models (Kerckhoffs) publish the active stress as
  // monitored value 0. Ionic models (ToRORdLand) do not: their Ta is a state
  // variable (49), so it has to be read directly.
  if (ta_var_index < 0)
    ephy.get_cells().get_monitored_values(0, ta);
  else
    ephy.get_cells().get_var(ta_var_index, ta);

  // Duas escalas, com papeis distintos:
  //   ta_scale (membro) = conversao kPa -> unidade de pressao do modelo
  //   ta_scale_global   = contratilidade global, -tascale
  // Ambas uniformes, aplicadas no campo NODAL. A terceira escala, a do
  // material, e por elemento e entra so na montagem.
  ta = ta * ((ta_scale != 0.0) ? ta_scale : T_ref) * ta_scale_global;

  // Tensao efetivamente entregue a montagem: a media nodal de Ta em cada
  // elemento vezes o ta_scale do material daquele elemento. Por elemento,
  // porque ta_scale e constante por elemento -- um no na fronteira entre
  // marcadores nao tem valor efetivo unico.
  elas.effective_active_tension(ta, ta_applied);

  update_active_stiffness();

  // Diagnostic only -- the values are NOT modified here. A negative active
  // tension is a compressive fibre stress, which is what inverts elements a
  // few steps later, so this line is the warning that precedes the negative
  // Jacobian.
  if (lambda_rate_on && ta.n_elem > 0 && ta.min() < 0.0)
    report_negative_ta();
}

void CardiacElectromechanic::update_active_stiffness()
{
  // The stabilisation only has a job when the active tension actually depends
  // on the shortening velocity. Without -lamrate, ZETAS and ZETAW stay at
  // rest, Ta has no velocity dependence, and there is no feedback loop to
  // stabilise -- so leave the mechanics exactly as it was.
  if (!stabilize_active || !lambda_rate_on || land_xs_index < 0
      || land_xw_index < 0 || lambda_raw.n_elem != ta.n_elem)
  {
    ka.reset();
    return;
  }

  if (land_xs.n_elem != ta.n_elem)
  {
    land_xs.zeros(ta.n_elem);
    land_xw.zeros(ta.n_elem);
    land_zetas.zeros(ta.n_elem);
    land_zetaw.zeros(ta.n_elem);
  }

  ephy.get_cells().get_var(land_xs_index, land_xs);
  ephy.get_cells().get_var(land_xw_index, land_xw);

  ka.set_size(ta.n_elem);
  for (arma::uword i = 0; i < ta.n_elem; i++)
    ka(i) = land_active_stiffness(lambda_raw(i), land_xs(i), land_xw(i));

  // Ka comes out of the Land submodel in kPa, exactly like Ta. It therefore
  // needs the SAME conversion that ta received two lines above; using a
  // different scale would leave a factor of 5000 between the stabilisation
  // term and the tension it is supposed to stabilise.
  ka = ka * ((ta_scale != 0.0) ? ta_scale : T_ref);
}

void CardiacElectromechanic::report_negative_ta()
{
  const int n_neg = static_cast<int>(arma::accu(ta < 0.0));

  cout << "  [MEF  ] *** WARNING: negative active tension at " << n_neg
       << "/" << ta.n_elem << " node(s), min = " << ta.min() << endl;

  // Without the Land state indices there is nothing further to say; the
  // decomposition below is specific to that submodel.
  if (land_xs_index < 0 || land_xw_index < 0
      || land_zetas_index < 0 || land_zetaw_index < 0)
    return;

  // Cells::get_var writes into v(i) without resizing it, exactly like the
  // pre-sized `ta` and `vm_node` buffers filled elsewhere. Sizing these here
  // rather than in config() keeps the allocation next to its only use and
  // costs nothing after the first call.
  if (land_xs.n_elem != ta.n_elem)
  {
    land_xs.zeros(ta.n_elem);
    land_xw.zeros(ta.n_elem);
    land_zetas.zeros(ta.n_elem);
    land_zetaw.zeros(ta.n_elem);
  }

  ephy.get_cells().get_var(land_xs_index,    land_xs);
  ephy.get_cells().get_var(land_xw_index,    land_xw);
  ephy.get_cells().get_var(land_zetas_index, land_zetas);
  ephy.get_cells().get_var(land_zetaw_index, land_zetaw);

  // Ta = Lfac*(Tref/dr)*[ (ZETAS+1)*XS + ZETAW*XW ]. The prefactor is
  // positive, so the sign is decided entirely by the bracket. Attribute each
  // negative node to whichever of the two terms is responsible:
  //
  //   "ZETAS": (ZETAS+1) < 0, i.e. the shortening was fast enough and
  //            sustained enough (tau_s = 24.9 ms) for the strong-crossbridge
  //            distortion to overshoot -1. This is the mechanism the old
  //            message assumed, and it needs |rate| > 1/249.1 = 0.004.
  //
  //   "ZETAW": (ZETAS+1) >= 0 but ZETAW*XW still outweighs (ZETAS+1)*XS.
  //            tau_w = 2.46 ms, so ZETAW tracks the instantaneous rate
  //            almost algebraically at any usable dtc, and when XS is small
  //            -- diastole, early systole, or any node not yet activated --
  //            a rate well inside -lamratemax is already enough. This is the
  //            common case, and the reason tightening -lamratemax alone does
  //            not fix the problem.
  int n_zetas = 0, n_zetaw = 0;
  double worst_rate = 0.0;
  double xs_at_worst = 0.0, xw_at_worst = 0.0;
  double worst_ta = 0.0;

  const bool have_rate = (lambda_rate_node.n_elem == ta.n_elem);

  for (arma::uword i = 0; i < ta.n_elem; i++)
  {
    if (ta(i) >= 0.0) continue;

    if (land_zetas(i) + 1.0 < 0.0) n_zetas++;
    else                           n_zetaw++;

    if (ta(i) < worst_ta)
    {
      worst_ta    = ta(i);
      xs_at_worst = land_xs(i);
      xw_at_worst = land_xw(i);
      worst_rate  = have_rate ? lambda_rate_node(i) : 0.0;
    }
  }

  cout << "  [MEF  ]   dominant term: ZETAW*XW at " << n_zetaw
       << " node(s), (ZETAS+1)<0 at " << n_zetas << " node(s)" << endl;
  cout << "  [MEF  ]   at the worst node: XS = " << xs_at_worst
       << ", XW = " << xw_at_worst
       << ", ZETAS = " << land_zetas.min()
       << ", ZETAW = " << land_zetaw.min() << endl;

  if (have_rate)
  {
    // The rate at which THIS node would have changed sign, from
    //     rate_crit = -XS / (249.1*XS + 24.6*XW).
    // Printing it next to the actual rate shows directly whether
    // -lamratemax could ever have prevented this node, or whether the
    // threshold was already below any cap worth using.
    const double denom = 249.128 * xs_at_worst + 24.639 * xw_at_worst;
    const double rate_crit = (denom > 0.0) ? -xs_at_worst / denom : 0.0;
    cout << "  [MEF  ]   d(lambda)/dt at that node = " << worst_rate
         << ", sign-change threshold = " << rate_crit << endl;

    if (n_zetaw > n_zetas && lambda_rate_max_clip > 0.0
        && std::fabs(rate_crit) < lambda_rate_max_clip)
      cout << "  [MEF  ]   the threshold is BELOW -lamratemax ("
           << lambda_rate_max_clip << "): tightening the clip alone"
           << " will not fix this." << endl;
  }
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

void CardiacElectromechanic::save_node_fields(int frame)
{
  // Both quantities live on the nodes: Cells is built with one ODE system
  // per mesh point (Eikonal::init), and the active tension is a nodal field
  // in the material (allocate_Ta(npoints)).

  // Membrane potential -- ionic models only. Kerckhoffs has no Vm, and its
  // state 0 is the contractile element length, which must not be written
  // into a field named vm.
  if (vm_var_index >= 0)
  {
    ephy.get_cells().get_var(vm_var_index, vm_node);
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

void CardiacElectromechanic::update_pv_jacobian(double t)
{
  const double dp = 0.05 * mmHg_to_p;   // 0.05 mmHg perturbation

  // baseline
  Solve_System(t, p_lv_cur, p_rv_cur);
  circ_solves++;
  double v0_lv = elas.volume_LV();
  double v0_rv = elas.volume_RV();

  // perturb LV pressure
  Solve_System(t, p_lv_cur + dp, p_rv_cur);
  circ_solves++;
  Jpv(0,0) = (elas.volume_LV() - v0_lv) / dp;
  Jpv(1,0) = (elas.volume_RV() - v0_rv) / dp;

  // perturb RV pressure
  Solve_System(t, p_lv_cur, p_rv_cur + dp);
  circ_solves++;
  Jpv(0,1) = (elas.volume_LV() - v0_lv) / dp;
  Jpv(1,1) = (elas.volume_RV() - v0_rv) / dp;

  Jpv_valid = true;

  if (!quiet_solve)
    cout << "    [PRESS] dV/dp = [" << Jpv(0,0) << " " << Jpv(0,1) << "; "
         << Jpv(1,0) << " " << Jpv(1,1) << "] mL per pressure unit" << endl;

  if (std::fabs(arma::det(Jpv)) < 1e-12)
  {
    cout << "   [coupling] WARNING: dV/dp is near singular; the two cavities"
         << " respond almost identically to both pressures." << endl;
  }
}

void CardiacElectromechanic::pressures_for_volumes(double Vlv_target,
                                                    double Vrv_target,
                                                    double t,
                                                    double & plv_mmHg,
                                                    double & prv_mmHg)
{
  // INVARIANT: this routine must NOT advance the cell model. It performs
  // several 3D solves per circulation step, and the 0D solver may call it
  // more than once per step; advancing ephy here would make the active
  // tension run ahead of the circulation clock by a factor equal to the
  // number of solves. `ta` is set once per step in solve_circulation() and
  // is held fixed throughout the inversion.
  const double tol_vol = 0.05;   // mL
  const int    max_it  = 8;

  int solves0 = circ_solves;

  if (!Jpv_valid)
  {
    if (!quiet_solve)
      cout << "    [PRESS] rebuilding dV/dp Jacobian (3 solves)" << endl;
    update_pv_jacobian(t);
  }

  arma::vec2 F;
  int it = 0;

  // evaluate at the current pressures
  Solve_System(t, p_lv_cur, p_rv_cur);
  circ_solves++;
  F(0) = elas.volume_LV() - Vlv_target;
  F(1) = elas.volume_RV() - Vrv_target;

  double f0 = arma::norm(F, 2);

  while (arma::norm(F, 2) > tol_vol && it < max_it)
  {
    arma::vec2 dp;
    bool ok = arma::solve(dp, Jpv, -F);
    if (!ok)
    {
      cout << "   [coupling] Jacobian solve failed; rebuilding" << endl;
      update_pv_jacobian(t);
      ok = arma::solve(dp, Jpv, -F);
      if (!ok) break;
    }

    // damped update: cavity pressures should not jump wildly in one step
    double maxstep = 20.0 * mmHg_to_p;
    double n = arma::norm(dp, 2);
    if (n > maxstep) dp *= maxstep / n;

    p_lv_cur += dp(0);
    p_rv_cur += dp(1);

    // Physiological bracket. If the 0D model asks for a volume the 3D model
    // cannot reach (too stiff, or a target outside its operating range),
    // Newton would otherwise run away to absurd pressures and take the whole
    // simulation with it.
    const double p_min = -20.0 * mmHg_to_p;
    const double p_max = 400.0 * mmHg_to_p;
    bool clipped = false;
    if (p_lv_cur < p_min) { p_lv_cur = p_min; clipped = true; }
    if (p_lv_cur > p_max) { p_lv_cur = p_max; clipped = true; }
    if (p_rv_cur < p_min) { p_rv_cur = p_min; clipped = true; }
    if (p_rv_cur > p_max) { p_rv_cur = p_max; clipped = true; }
    if (clipped)
    {
      cout << "   [coupling] WARNING t=" << t << ": pressure hit the"
           << " physiological bracket. The requested volumes (LV "
           << Vlv_target << " mL, RV " << Vrv_target
           << " mL) may be unreachable for this 3D model." << endl;
      Jpv_valid = false;
      break;
    }

    Solve_System(t, p_lv_cur, p_rv_cur);
    circ_solves++;
    F(0) = elas.volume_LV() - Vlv_target;
    F(1) = elas.volume_RV() - Vrv_target;

    it++;
  }

  // If Newton stalled, the Jacobian is stale: force a rebuild next time.
  if (arma::norm(F, 2) > tol_vol)
  {
    Jpv_valid = false;
    cout << "   [coupling] t=" << t << " residual |dV| = "
         << arma::norm(F, 2) << " mL after " << it
         << " iterations (target " << tol_vol << "); Jacobian will be rebuilt"
         << endl;
  }
  else if (it > 3 || arma::norm(F,2) > 0.5 * f0)
  {
    Jpv_valid = false;   // slow convergence -> refresh
  }

  last_newton_its = it;
  last_dV         = arma::norm(F, 2);
  last_solves     = circ_solves - solves0;

  plv_mmHg = p_lv_cur / mmHg_to_p;
  prv_mmHg = p_rv_cur / mmHg_to_p;
}

void CardiacElectromechanic::solve_circulation(int num_beats, double dt_circ)
{
  cout << "\n=== Closed-loop 3D-0D circulation (Regazzoni 2022) ===" << endl;

  timer.enter("solve_circulation");

  // config() ja resolveu -beats contra o teto imposto por -t. O argumento
  // chega direto da linha de comando pelo app e pode estar desatualizado.
  if (circ_num_beats > 0 && num_beats != circ_num_beats)
  {
    cout << " *** -beats " << num_beats << " ajustado para " << circ_num_beats
         << " (limite imposto por -t)" << endl;
    num_beats = circ_num_beats;
  }

  // --- initialisation, mirroring solve() ------------------------------
  // pre_solve() builds the elasticity matrices/vectors and the reference
  // configuration; initial_conditions() resets the cell model state.
  // Skipping either leaves the linear system unassembled, which shows up as
  // a preconditioner failure (PETSc reason -11) on the very first solve.
  cout << " Time step of mechanics: " << dt_mech << endl;
  elas.pre_solve();
  ephy.initial_conditions();

  // active tension at t = 0
  update_active_tension();

  // --- initialise the 0D model from the actual 3D cavity volumes ------
  // Sobrescrever V_LV/V_RV com os volumes da malha MUDA o volume sanguineo
  // estressado: o estado inicial publicado tem V_LV = 118.5 e V_RV = 166.2
  // mL, e uma malha tipica esta dezenas de mL abaixo disso. Sem compensacao,
  // o modelo inteiro roda subenchido e o efeito e invisivel no log.
  const double vtot0 = CommandLineArgs::read("-circ_vtot",
      circ.total_volume(),
      "volume sanguineo estressado em mL (default: o do estado inicial)",
      "circulacao");

  double v_lv0 = elas.volume_LV();
  double v_rv0 = elas.volume_RV();
  circ.set_initial_volumes(v_lv0, v_rv0);
  const double dv_realloc = circ.rebalance_total_volume(vtot0);

  cout << " Initial cavity volumes from the 3D mesh: LV = " << v_lv0
       << " mL, RV = " << v_rv0 << " mL" << endl;
  cout << " Stressed blood volume held at " << vtot0 << " mL ("
       << dv_realloc << " mL reallocated to the venous reservoirs)" << endl;
  cout << " Pressure conversion: 1 mmHg = " << mmHg_to_p
       << " model pressure units" << endl;

  // --- register the 3D coupling --------------------------------------
  circ.set_BiV_coupling(
    [this](double Vlv, double Vrv, double t, double & plv, double & prv)
    {
      this->pressures_for_volumes(Vlv, Vrv, t, plv, prv);
    });

  // --- output ---------------------------------------------------------
  circ_file.open((filename + string("_circulation.dat")).c_str());
  if (circ_file.is_open())
    circ_file << Regazzoni2020::history_header() << "\n";

  const double T = circ.period();
  int    nsteps  = static_cast<int>(std::round(T / dt_circ));
  double t = 0.0;

  // ---- cell-model sub-cycling ----------------------------------------
  // The circulation model works in SECONDS. The electrophysiology may be
  // configured in seconds or in milliseconds depending on the input file,
  // so the unit is DETECTED rather than assumed: the ephy final time should
  // span a whole number of beats, which pins the scale.
  //
  // The cell model is advanced ONCE per circulation step, never inside the
  // pressure inversion: pressures_for_volumes() runs several 3D solves per
  // step and advancing there would make Ta race ahead of the circulation by
  // a factor equal to the number of Newton iterations.
  double dt_ephy   = ephy.get_time_parameters().get_dt();
  double stop_ephy = ephy.get_time_parameters().get_stop();

  if (ephy_time_per_second <= 0.0)
  {
    // auto-detect: is the ephy span closer to T seconds or to T*1000 ms?
    double err_s  = std::fabs(stop_ephy - T);
    double err_ms = std::fabs(stop_ephy - T * 1000.0);
    ephy_time_per_second = (err_s <= err_ms) ? 1.0 : 1000.0;

    cout << " Cell model span = " << stop_ephy
         << ", heart period = " << T << " s -> electrophysiology is in "
         << (ephy_time_per_second == 1.0 ? "SECONDS" : "MILLISECONDS")
         << endl;

    double rel = std::min(err_s / T, err_ms / (T * 1000.0));
    if (rel > 0.05)
      cout << " *** WARNING: the cell-model span does not match a whole"
           << " number of beats; unit detection is unreliable. Set it"
           << " explicitly with set_ephy_time_per_second()." << endl;
  }

  double dt_circ_ephy = dt_circ * ephy_time_per_second;
  int n_sub = static_cast<int>(std::round(dt_circ_ephy / dt_ephy));
  if (n_sub < 1) n_sub = 1;

  double mismatch = std::fabs(n_sub * dt_ephy - dt_circ_ephy);
  cout << " Cell model: dt = " << dt_ephy << ", circulation step = "
       << dt_circ_ephy << " (same units) -> " << n_sub
       << " sub-step(s) per step" << endl;
  if (mismatch > 1e-9 * dt_circ_ephy)
  {
    cout << " *** WARNING: the circulation step is not a multiple of the"
         << " cell-model step." << endl;
    cout << "     " << n_sub << " x " << dt_ephy << " = "
         << n_sub * dt_ephy << " vs " << dt_circ_ephy
         << "  (drift of " << mismatch << " per step)." << endl;
  }

  cout << " Heart period = " << T << " s, dt = " << dt_circ
       << " s -> " << nsteps << " steps/beat, " << num_beats << " beat(s)"
       << endl;

  // Period of the LAT stimulus, in the ELECTROPHYSIOLOGY time unit: lat,
  // stim_duration and the cell-model clock all live in that unit, which may
  // be ms or s (see ephy_time_per_second above), while T is in seconds.
  // Without this the stimulus window would be tested against the wrong scale.
  const double stim_period = T * ephy_time_per_second;

  // ---- o relogio da EP tem de cobrir a rodada INTEIRA -----------------
  // Se estiver curto, estica: no caminho acoplado quem manda no tempo e este
  // laco, nao o -t. A excecao abaixo so dispara se a aritmetica ficar
  // inconsistente, e e MUITO preferivel ao no-op silencioso de
  // Eikonal::advance(), que congela as celulas sem emitir nada.
  const double total_ephy = num_beats * stim_period;
  if (stop_ephy < total_ephy - 0.5 * dt_ephy)
  {
    const double novo = total_ephy + 2.0 * dt_ephy;
    cout << " EP span esticado de " << stop_ephy << " para " << novo
         << " para cobrir " << num_beats << " batimento(s)" << endl;
    ephy.time_parameters().set_stop(novo);
    stop_ephy = ephy.get_time_parameters().get_stop();
  }
  if (stop_ephy < total_ephy - 0.5 * dt_ephy)
    throw std::runtime_error("EP span shorter than the run: the cell model "
                             "would silently stop integrating partway "
                             "through. Check -t, -hr and -beats.");

  if (stim_amplitude != 0.0)
    cout << " LAT stimulus: amplitude = " << stim_amplitude
         << ", duration = " << stim_duration
         << ", repeating every " << stim_period
         << " (cell-model time units)" << endl;

  // The data writer was sized and opened in config(); here we only use it.
  int out_every = circ_out_every;
  int n_frames  = circ_n_frames;

  int frame = 0;
  elas.output_vtk(0, frame);
  save_node_fields(frame);
  cout << " Output: saving every " << out_every << " step(s), "
       << n_frames << " frame(s) allocated" << endl;
  cout << " Log: detailed report every " << circ_log_every << " step(s)"
       << endl;
  cout << "\n Legend: [ODE] cell model | [CIRC] 0D circulation |"
       << " [PRESS] pressure inversion | [STATE] volumes/pressures |"
       << " [FLOW] valve flows" << endl;

  for (int beat = 0; beat < num_beats; beat++)
  {
   cout << "\n--- beat " << beat + 1 << "/" << num_beats << " ---" << endl;

   for (int i = 0; i < nsteps; i++)
   {
    bool log_now = (i % circ_log_every == 0) || (i == nsteps - 1);
    quiet_solve = !log_now;

    if (log_now)
    {
      // std::defaultfloat restores the FORMAT but not the PRECISION, so the
      // setprecision(1) below used to leak into every line printed until
      // something else set it again -- which is why the MEF diagnostics came
      // out as "lambda_f: min = 1, mean = 1, max = 1" and hid exactly the
      // digits they exist to show. Save and restore it explicitly.
      const std::streamsize prec0 = cout.precision();
      cout << "\n[t = " << std::fixed << std::setprecision(1)
           << (t + dt_circ) * 1000.0 << " ms | step " << i + 1 << "/" << nsteps
           << " | beat " << beat + 1 << "/" << num_beats << "]"
           << std::defaultfloat << endl;
      cout.precision(prec0);
    }

    // ---- 0. mechano-electric feedback -----------------------------
    // lambda_f comes from the deformation gradient left by the LAST
    // converged mechanical solve of the previous step, so it must be read
    // before the cells advance. It is read once per circulation step, never
    // inside pressures_for_volumes(): that routine runs several trial
    // solves, and the stretch of a rejected trial pressure has no meaning.
    update_fiber_stretch(dt_circ_ephy);

    // ---- 1. cell model -------------------------------------------
    // advance over one circulation step, then freeze: ta stays constant
    // while the pressure inversion iterates.
    // "stimulus on" must reflect whether any node was actually stimulated
    // in this step, not merely that stimulation is enabled.
    bool stimulating = false;
    const arma::vec & lat_v = ephy.get_lat();
    for (int k = 0; k < n_sub; k++)
    {
      if (stim_amplitude != 0.0)
      {
        double tk = ephy.get_time_parameters().time();
        if (stim_period > 0.0) tk = std::fmod(tk, stim_period);
        if (arma::any(tk >= lat_v && tk < lat_v + stim_duration))
          stimulating = true;
        ephy.apply_lat_stimulus(stim_amplitude, stim_duration, stim_period);
      }
      ephy.advance();
    }
    update_active_tension();

    // A stiff ionic model integrated with too large a step blows up to NaN
    // in the membrane potential first; Ta often keeps decaying smoothly for
    // a while, so the run looks alive while the electrophysiology is dead.
    if (!ta.is_finite())
    {
      cout << "\n *** the active tension contains NaN/Inf at t = "
           << (t + dt_circ) * 1000.0 << " ms." << endl;
      cout << "     The cell model has diverged. With explicit Euler the"
           << " ToRORd model needs dt <= ~0.002 ms;" << endl;
      cout << "     reduce -dt (0.001 is safe) or use an implicit solver."
           << endl;
      break;
    }

    if (log_now)
    {
      // Duas grandezas diferentes, e a distincao importa: 'ta' e o valor
      // BRUTO do modelo celular e nao muda quando se mexe no ta_scale do
      // XML; 'ta_applied' e o que a montagem realmente usa.
      // 'nodal'   = ja inclui kPa->unidades e o -tascale global, uniforme.
      // 'applied' = tambem inclui o ta_scale do material, por elemento.
      // Nao sao proporcionais: o maximo nodal pode cair num no cujos
      // elementos vizinhos sejam todos passivos.
      cout << "  [ODE  ] " << n_sub << " sub-steps"
           << (stimulating ? " (stimulus on)" : "")
           << " -> Ta max = " << ta.max() << " nodal";
      if (ta_applied.n_elem > 0)
        cout << ", " << ta_applied.max() << " applied";
      cout << endl;
    }

    if (i == 0 && beat == 0 && arma::max(arma::abs(ta)) == 0.0)
    {
      cout << " *** WARNING: the active tension is zero at the start of"
           << " the beat." << endl;
      cout << "     With Kerckhoffs, Ta is non-zero only for"
           << " lat <= t <= lat + 0.483 s, where lat comes from the"
           << " per-node LAT in the mesh or, if absent, from the"
           << " <pvloop passive_time=...> attribute." << endl;
      cout << "     Check that passive_time is smaller than the heart"
           << " period (" << T << " s), otherwise the ventricle will only"
           << " fill passively and never eject." << endl;
    }

    // ---- 2. circulation + 3. pressure inversion -------------------
    // circ.step() calls pressures_for_volumes() internally, which runs the
    // 3D solves needed to match the volumes the 0D model asks for.
    if (log_now)
      cout << "  [CIRC ] advancing 0D model, dt = " << dt_circ << " s" << endl;

    circ.step(t, dt_circ);
    t += dt_circ;

    if (log_now)
    {
      cout << "  [PRESS] " << last_newton_its << " Newton it, |dV| = "
           << last_dV << " mL, " << last_solves << " 3D solve(s)" << endl;
      cout << "  [STATE] LV: V = " << circ[Regazzoni2020::V_LV]
           << " mL, p = " << circ.pressure_LV() << " mmHg"
           << "  |  RV: V = " << circ[Regazzoni2020::V_RV]
           << " mL, p = " << circ.pressure_RV() << " mmHg" << endl;
      cout << "  [FLOW ] MV = " << circ.flow_MV() << ", AV = " << circ.flow_AV()
           << ", TV = " << circ.flow_TV() << ", PV = " << circ.flow_PV()
           << " mL/s  | Vtot = " << circ.total_volume() << " mL" << endl;
    }

    if (circ_file.is_open())
    {
      circ_file << circ.history_row(t) << "\n";
      circ_file.flush();
    }

    if (i % out_every == 0 && frame < n_frames)
    {
      frame++;
      elas.output_vtk(0, frame);
      save_node_fields(frame);
    }
   }
  }

  quiet_solve = false;

  if (circ_file.is_open())
    circ_file.close();

  cout << "\n Total 3D solves: " << circ_solves
       << " (" << double(circ_solves) / (num_beats * nsteps)
       << " per circulation step)" << endl;
  cout << " Final blood volume: " << circ.total_volume() << " mL" << endl;

  timer.leave();
}

// ======================================================================
//  PV loop POR FASES  (Wang et al. 2025, eLife 14:RP106555)
// ======================================================================
//
//  A diferenca estrutural em relacao a solve_circulation() esta em UMA
//  linha do laco: aqui a pressao do passo sai de uma expressao explicita
//  (PiecewisePVLoop::advance_pressure) e vai direto para a mecanica, em
//  vez de ser a incognita de um Newton 2x2 sobre o volume.
//
//  Consequencias praticas:
//
//    - custo: UMA solucao 3D por passo, contra 4-8 do caminho -circ;
//    - estabilidade: o esquema e explicito no acoplamento p <-> V, entao
//      as fases isovolumetricas sao mantidas por penalidade e nao por
//      restricao. Ganho alto demais nas penalidades oscila; ganho baixo
//      demais deixa o volume derivar. Ver as notas em piecewise_pvloop.hpp;
//    - circulacao: ABERTA. Um Windkessel de dois elementos por ventriculo,
//      sem retorno venoso e sem atrios. Nao ha volume sanguineo total a
//      conservar, o que elimina a questao do volume estressado que domina
//      o transiente do caminho -circ, mas tambem elimina o mecanismo de
//      Frank-Starling entre batimentos.
//
//  Ordem dentro do passo, e por que ela e essa:
//
//    0. lambda_f vem do gradiente de deformacao deixado pela ULTIMA
//       solucao convergida, logo tem de ser lido ANTES das celulas;
//    1. celulas avancam n_sub sub-passos e Ta e congelada;
//    2. o controlador calcula p_LV/p_RV a partir de V e dV/dt do passo
//       anterior;
//    3. UMA solucao mecanica com essas pressoes e essa Ta;
//    4. os volumes resultantes voltam para o controlador, que atualiza
//       dV/dt e avalia as transicoes de fase.

void CardiacElectromechanic::solve_piecewise(int num_beats, double dt_pw_ms)
{
  cout << "\n=== PV loop por fases (Wang et al. 2025) ===" << endl;

  timer.enter("solve_piecewise");

  // config() ja resolveu -beats e -dtp contra o teto imposto por -t. Os
  // argumentos chegam direto da linha de comando pelo app e podem estar
  // desatualizados; ESTES membros sao a verdade.
  if (pw_num_beats > 0 && num_beats != pw_num_beats)
  {
    cout << " *** -beats " << num_beats << " ajustado para " << pw_num_beats
         << " (limite imposto por -t)" << endl;
    num_beats = pw_num_beats;
  }
  if (std::fabs(dt_pw_ms - pw_dt_ms) > 1e-12)
    dt_pw_ms = pw_dt_ms;

  // --- inicializacao, espelhando solve_circulation() -------------------
  // pre_solve() monta as matrizes e a configuracao de referencia;
  // initial_conditions() reinicia o estado celular. Pular qualquer um dos
  // dois deixa o sistema linear sem montagem, o que aparece como falha de
  // precondicionador (PETSc reason -11) ja na primeira solucao.
  cout << " Time step of mechanics: " << dt_mech << endl;
  elas.pre_solve();
  ephy.initial_conditions();

  update_active_tension();

  // --- estado inicial do controlador -----------------------------------
  const double v_lv0 = elas.volume_LV();
  const double v_rv0 = elas.volume_RV();
  pw.init(v_lv0, v_rv0);

  cout << " Volumes iniciais da malha: VE = " << v_lv0
       << " mL, VD = " << v_rv0 << " mL" << endl;
  cout << " Conversao de pressao: 1 kPa = " << kPa_to_p
       << " unidades do modelo" << endl;

  if (v_rv0 > v_lv0)
    cout << " *** NOTA: o volume da cavidade DIREITA excede o da esquerda"
         << " nesta malha. Confirme os marcadores (30 = VE, 20 = VD) antes"
         << " de calibrar os limiares de ejecao." << endl;

  // --- saida -----------------------------------------------------------
  pw_file.open((filename + string("_piecewise.dat")).c_str());
  if (pw_file.is_open())
    pw_file << PiecewisePVLoop::history_header() << "\n";

  // --- relogios --------------------------------------------------------
  //
  // Tres unidades de tempo convivem neste laco e a maior parte dos erros
  // caros mora exatamente aqui:
  //
  //   ms                : o controlador e a Tabela 1 do artigo;
  //   unidade do solver : o relogio da EP, s ou ms conforme -tunit;
  //   unidade do modelo : a do proprio modelo celular, tratada por Cells.
  //
  // ephy_time_per_second foi fixado em config() a partir de -tunit, entao
  // aqui ele NAO e autodetectado -- ao contrario do caminho -circ, que
  // aceita malhas antigas sem essa declaracao.
  const double units_per_ms = (ephy_time_per_second > 0.0)
                            ? ephy_time_per_second / 1000.0 : 1.0;

  const double dt_ephy   = ephy.get_time_parameters().get_dt();
  const double dt_pw_ep  = dt_pw_ms * units_per_ms;   // passo em unid. do solver

  int n_sub = static_cast<int>(std::round(dt_pw_ep / dt_ephy));
  if (n_sub < 1) n_sub = 1;

  const double mismatch = std::fabs(n_sub * dt_ephy - dt_pw_ep);
  cout << " Cell model: dt = " << dt_ephy << ", passo do controlador = "
       << dt_pw_ep << " (mesmas unidades) -> " << n_sub
       << " sub-passo(s) por passo" << endl;
  if (mismatch > 1e-9 * dt_pw_ep)
  {
    cout << " *** ATENCAO: o passo do controlador nao e multiplo do passo"
         << " do modelo celular." << endl;
    cout << "     " << n_sub << " x " << dt_ephy << " = " << n_sub * dt_ephy
         << " vs " << dt_pw_ep << " (deriva de " << mismatch
         << " por passo)." << endl;
  }

  const double RR_ms      = pw.timing().RR;
  const double total_ms   = pw.total_time(num_beats);
  const int    nsteps     = static_cast<int>(std::round(total_ms / dt_pw_ms));

  // Estimulo LAT: periodo e deslocamento na unidade do SOLVER, porque e
  // contra o relogio da EP que Eikonal::apply_lat_stimulus os compara.
  // O deslocamento e o fim da fase 1: no artigo, "Phase 2: Isovolumetric
  // contraction. Electrical activation ensues". Sem ele, a ativacao cairia
  // em t = 0, isto e, no meio da rampa de inflacao passiva.
  const double stim_period_ep = RR_ms * units_per_ms;
  const double stim_offset_ep = pw.stimulus_offset() * units_per_ms;

  cout << " Rodada: " << total_ms << " ms = " << nsteps << " passos de "
       << dt_pw_ms << " ms (" << num_beats << " batimento(s) de " << RR_ms
       << " ms apos a fase 0)" << endl;

  if (stim_amplitude != 0.0)
    cout << " Estimulo LAT: amplitude = " << stim_amplitude
         << ", duracao = " << stim_duration
         << ", periodo = " << stim_period_ep
         << ", offset = " << stim_offset_ep
         << " (unidades do solver)" << endl;
  else
    cout << " Sem estimulo LAT: o modelo celular e dirigido diretamente"
         << " pelo LAT (Kerckhoffs). A ativacao NAO respeita o offset de"
         << " fase; considere ajustar o LAT da malha em "
         << pw.stimulus_offset() << " ms." << endl;

  // --- o relogio da EP tem de cobrir a rodada inteira -------------------
  double stop_ephy = ephy.get_time_parameters().get_stop();
  const double total_ep = total_ms * units_per_ms;
  if (stop_ephy < total_ep - 0.5 * dt_ephy)
  {
    const double novo = total_ep + 2.0 * dt_ephy;
    cout << " EP span esticado de " << stop_ephy << " para " << novo
         << " para cobrir a rodada" << endl;
    ephy.time_parameters().set_stop(novo);
    stop_ephy = ephy.get_time_parameters().get_stop();
  }
  if (stop_ephy < total_ep - 0.5 * dt_ephy)
    throw std::runtime_error("EP span menor que a rodada: o modelo celular "
                             "pararia de integrar no meio, sem avisar. "
                             "Verifique -t, -hr/-pwrr, -pwinit e -beats.");

  // --- primeiro frame ---------------------------------------------------
  int frame = 0;
  elas.output_vtk(0, frame);
  save_node_fields(frame);
  cout << " Output: salvando a cada " << pw_out_every << " passo(s), "
       << pw_n_frames << " frame(s) alocados" << endl;
  cout << " Log: relatorio detalhado a cada " << pw_log_every << " passo(s)"
       << endl;
  cout << "\n Legenda: [ODE] modelo celular | [PW] controlador de fases |"
       << " [MECH] mecanica 3D" << endl;

  double t_ms   = 0.0;
  int    solves = 0;
  int    beat_reported = -1;

  for (int i = 0; i < nsteps; i++)
  {
    const bool log_now = (i % pw_log_every == 0) || (i == nsteps - 1);
    quiet_solve = !log_now;

    const double t_end_ms = t_ms + dt_pw_ms;

    const int beat = pw.beat_of(t_end_ms);
    if (beat != beat_reported && beat >= 0)
    {
      cout << "\n--- batimento " << beat + 1 << "/" << num_beats << " ---"
           << endl;
      beat_reported = beat;
    }

    if (log_now)
    {
      const std::streamsize prec0 = cout.precision();
      cout << "\n[t = " << std::fixed << std::setprecision(1) << t_end_ms
           << " ms | passo " << i + 1 << "/" << nsteps
           << " | VE " << PiecewisePVLoop::phase_name(pw.state_LV().phase)
           << " | VD " << PiecewisePVLoop::phase_name(pw.state_RV().phase)
           << "]" << std::defaultfloat << endl;
      cout.precision(prec0);
    }

    // ---- 0. realimentacao mecano-eletrica ---------------------------
    update_fiber_stretch(dt_pw_ep);

    // ---- 1. modelo celular ------------------------------------------
    bool stimulating = false;
    const arma::vec & lat_v = ephy.get_lat();
    for (int k = 0; k < n_sub; k++)
    {
      if (stim_amplitude != 0.0)
      {
        double tk = ephy.get_time_parameters().time() - stim_offset_ep;
        if (stim_period_ep > 0.0)
          tk = std::fmod(std::fmod(tk, stim_period_ep) + stim_period_ep,
                         stim_period_ep);
        if (tk >= 0.0 && arma::any(tk >= lat_v && tk < lat_v + stim_duration))
          stimulating = true;

        ephy.apply_lat_stimulus(stim_amplitude, stim_duration,
                                stim_period_ep, stim_offset_ep);
      }
      ephy.advance();
    }
    update_active_tension();

    if (!ta.is_finite())
    {
      cout << "\n *** a tensao ativa contem NaN/Inf em t = " << t_end_ms
           << " ms." << endl;
      cout << "     O modelo celular divergiu. Com Euler explicito o ToRORd"
           << " precisa de dt <= ~0.002 ms;" << endl;
      cout << "     reduza -dt ou use um solver implicito." << endl;
      break;
    }

    if (log_now)
    {
      cout << "  [ODE  ] " << n_sub << " sub-passos"
           << (stimulating ? " (estimulo ligado)" : "")
           << " -> Ta max = " << ta.max() << " nodal";
      if (ta_applied.n_elem > 0)
        cout << ", " << ta_applied.max() << " applied";
      cout << endl;
    }

    // ---- 2. controlador de fases ------------------------------------
    pw.advance_pressure(t_end_ms, dt_pw_ms);

    const double p_lv_kPa = pw.pressure_LV();
    const double p_rv_kPa = pw.pressure_RV();

    if (log_now)
      cout << "  [PW   ] p_VE = " << p_lv_kPa << " kPa, p_VD = " << p_rv_kPa
           << " kPa  (P_art VE = " << pw.state_LV().P_art << " kPa)" << endl;

    // ---- 3. UMA solucao mecanica -------------------------------------
    // A conversao kPa -> unidade de pressao do modelo e a MESMA que a
    // tensao ativa recebe (kPa_to_p), de proposito: um fator de 1000
    // entre a pressao da cavidade e a tensao ativa seria quase invisivel
    // nos resultados e fatal na fisiologia.
    p_lv_cur = p_lv_kPa * kPa_to_p;
    p_rv_cur = p_rv_kPa * kPa_to_p;

    Solve_System(t_end_ms, p_lv_cur, p_rv_cur);
    solves++;

    // ---- 4. volumes de volta ao controlador --------------------------
    const double v_lv = elas.volume_LV();
    const double v_rv = elas.volume_RV();

    pw.update_volumes(t_end_ms, dt_pw_ms, v_lv, v_rv);

    t_ms = t_end_ms;

    if (log_now)
      cout << "  [MECH ] V_VE = " << v_lv << " mL (dV/dt = "
           << pw.state_LV().Vdot << " mL/ms)"
           << " | V_VD = " << v_rv << " mL (dV/dt = "
           << pw.state_RV().Vdot << " mL/ms)" << endl;

    if (pw_file.is_open())
    {
      pw_file << pw.history_row(t_ms) << "\n";
      pw_file.flush();
    }

    if (i % pw_out_every == 0 && frame < pw_n_frames - 1)
    {
      frame++;
      elas.output_vtk(0, frame);
      save_node_fields(frame);
    }
  }

  quiet_solve = false;

  if (pw_file.is_open())
    pw_file.close();

  // --- resumo -----------------------------------------------------------
  cout << "\n Total de solucoes 3D: " << solves
       << " (" << (nsteps > 0 ? double(solves) / nsteps : 0.0)
       << " por passo)" << endl;

  const PiecewisePVLoop::ChamberState & slv = pw.state_LV();
  const PiecewisePVLoop::ChamberState & srv = pw.state_RV();

  cout << " VE: V_ED = " << slv.V_ED << " mL, V_ES = " << slv.V_ES
       << " mL, VS = " << pw.stroke_volume_LV() << " mL, FE = "
       << pw.ejection_fraction_LV() << " %" << endl;
  cout << " VD: V_ED = " << srv.V_ED << " mL, V_ES = " << srv.V_ES
       << " mL, VS = " << pw.stroke_volume_RV() << " mL, FE = "
       << pw.ejection_fraction_RV() << " %" << endl;

  // O teto de |dP| e uma rede de seguranca. Se ele atuou muito, o ganho
  // das penalidades esta alto demais para este dt e a fase isovolumetrica
  // esta sendo suavizada de um jeito que a calibracao nao ve.
  if (slv.n_clipped > 0 || srv.n_clipped > 0)
    cout << " *** ATENCAO: o teto -pwdpmax atuou " << slv.n_clipped
         << " vez(es) no VE e " << srv.n_clipped << " no VD. Reduza"
         << " -pwcpiso/-pwcviso ou -dtp antes de interpretar o laco."
         << endl;

  if (slv.V_ES >= slv.V_ED)
    cout << " *** ATENCAO: V_ES >= V_ED no VE: nao houve ejecao liquida."
         << " Verifique -tascale, o limiar -pwpej e se a ativacao caiu no"
         << " fim da fase 1." << endl;

  timer.leave();

  elas.timer.summary();
  timer.summary();
}