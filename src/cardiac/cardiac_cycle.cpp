#include <armadillo>
#include <cmath>
#include "util/command_line_args.h"
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
  double T_cli  = CommandLineArgs::read("-t", -1.0);
  double dt_cli = CommandLineArgs::read("-dt", -1.0);

  if (T_cli > 0.0 && T_cli != T)
  {
    cout << " total_time from the mesh (" << T << ") overridden by -t "
         << T_cli << endl;
    T = T_cli;
  }
  if (dt_cli > 0.0)
  {
    cout << " dt from the pvloop curve (" << dt_from_curve
         << ") overridden by -dt " << dt_cli << endl;
    dt = dt_cli;
  }
  else
    dt = dt_from_curve;

  dt_mech = dt;
  cout << "size: " << size << endl;
  cout << "total time: " << T << endl;
  cout << "Mechanical time step: " << dt_mech << endl;

  ephy.setup(monofile, cellmodel, odesolver, dt, T, pr, pa);
  ephy.update_matrix(false);
  ephy.init();

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
  if (T_ms < 50.0 || T_ms > 20000.0)
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
    // Land gives Ta in kPa; the x5 factor follows the fenicsx-pulse demo.
    set_active_tension_source(49, 5.0 * kPa_to_p);
    set_vm_source(0);                     // Vm is state variable 0
    double amp = CommandLineArgs::read("-stimamp", -53.0);
    double dur = CommandLineArgs::read("-stimdur", 2.0 * ephy.ms_to_solver_time());
    set_lat_stimulus(amp, dur);
    cout << " Active tension: state variable 49 (kPa) x 5.0 x " << kPa_to_p
         << " -> " << 5.0 * kPa_to_p << " model pressure units per kPa"
         << endl;
    cout << " LAT stimulus: amplitude " << amp << ", duration " << dur
         << " (solver time units)" << endl;
  }
  else
  {
    set_active_tension_source(-1, 0.0);   // monitored 0, scaled by T_ref
    set_vm_source(-1);                    // no membrane potential to save
    set_lat_stimulus(0.0, 0.002);         // no stimulus needed
  }


  elas.config(mshfile, parfile);
  elas.set_output_step(false);
  elas.init();

  // Size the data writer here, once. For the closed-loop circulation path
  // the number of frames follows from the run parameters; the length of the
  // pressure curve (curr_time) is irrelevant there because that curve is
  // not used at all.
  int use_circ  = CommandLineArgs::read("-circ", 0);
  circ_out_every = CommandLineArgs::read("-prc", 8);
  circ_log_every = CommandLineArgs::read("-logc", 25);
  if (circ_log_every < 1) circ_log_every = 1;
  if (circ_out_every < 1) circ_out_every = 1;

  if (use_circ)
  {
    int    num_beats = CommandLineArgs::read("-beats", 1);
    double dtc       = CommandLineArgs::read("-dtc", 1.0e-3);
    int nsteps = static_cast<int>(std::round(circ.period() / dtc));
    circ_n_frames = num_beats * nsteps / circ_out_every + 1;
    cout << " Circulation output: " << num_beats * nsteps
         << " steps, saving every " << circ_out_every << " -> "
         << circ_n_frames << " frame(s)" << endl;
    elas.setup_data_writer(circ_n_frames);
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

  //"Solve" the eikonal, to discover the lat in each element
  ephy.solve(basename);
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

  ta = ta * ((ta_scale != 0.0) ? ta_scale : T_ref);
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
  elas.store_point_field(frame, ta, string("active_stress"));
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
  double v_lv0 = elas.volume_LV();
  double v_rv0 = elas.volume_RV();
  circ.set_initial_volumes(v_lv0, v_rv0);

  cout << " Initial cavity volumes from the 3D mesh: LV = " << v_lv0
       << " mL, RV = " << v_rv0 << " mL" << endl;
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
      cout << "\n[t = " << std::fixed << std::setprecision(1)
           << (t + dt_circ) * 1000.0 << " ms | step " << i + 1 << "/" << nsteps
           << " | beat " << beat + 1 << "/" << num_beats << "]"
           << std::defaultfloat << endl;
    }

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
      cout << "  [ODE  ] " << n_sub << " sub-steps"
           << (stimulating ? " (stimulus on)" : "")
           << " -> Ta max = " << ta.max() << endl;

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