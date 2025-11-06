#include <armadillo>
#include "util/command_line_args.h"
#include "cardiac_cycle.hpp"
#include "util/pugixml.hpp"

CardiacElectromechanic::CardiacElectromechanic(const std::string &epmodel) : dt_mech(1.0), timer(), Ta0(0.0), dTa0(0.0)
{
  // TODO: change from EPHY to EP
  // create electrophysiology model
  // if (epmodel == "mono")
  //  ep = new MonodomainDeformation();
  // else
  //  ep = new BidomainDeformation();
}

CardiacElectromechanic::~CardiacElectromechanic()
{
  // freeing memory
  /*
  for (uint i = 0; i < vec_stress.size(); i++)
    delete vec_stress[i];
  for (uint i = 0; i < vec_fib.size(); i++)
    delete vec_fib[i];
  for (uint i = 0; i < vec_fib0.size(); i++)
    delete vec_fib0[i];

  for (uint i = 0; i < Ta_list.size(); i++)
    delete Ta_list[i];
  for (uint i = 0; i < p_lv.size(); i++)
    delete p_lv[i];
  for (uint i = 0; i < p_art.size(); i++)
    delete p_art[i];
  for (uint i = 0; i < p_ven.size(); i++)
    delete p_ven[i];
  for (uint i = 0; i < p_LA.size(); i++)
    delete p_LA[i];
  */
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
  cellmodel = CommandLineArgs::read("-c", "NP");
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

  //  cout <<"PV-loop parameters:\n";
  //  cout <<"C_art: " << C_art << " p_art: " << part << " R_per: " << R_per << " P_o: " << P_o << " Stroke volume: " << stroke_volume << "\n";

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
    // std::pair<double, double> par;
    // par.first=press;
    // par.second= (ta<1e-16) ? 0 : ta;
    // Ta_list.push_back(ta); //! Comentei aqui para testar
    p_lv.push_back(press);
    p_rv.push_back(press2);

    cout << "Pressure: " << press << " " << press2 << " Ta: " << ta << endl;
  }


  cout<<"-- Calculating active tension for all time steps --" <<endl; 
  lc = 1.9;

 /*  double max_ta = 0.0; 
  for (size_t i = 0; i < curr_time.size(); ++i)
  {
    double Ta_calculated = solveTa((curr_time[i] / 1000.0) - 0.136, dt / 1000.0);
    Ta_list.push_back(Ta_calculated);
    if(Ta_calculated > max_ta)
      max_ta = Ta_calculated; 
  }

  //Normalizando e multiplicando pelo _tref
  if (max_ta != 0.0) {
    for (double& ta : Ta_list) {
        ta = (ta * T_ref) / max_ta;
    }
  } else {
      std::cerr << "Warning: max_ta is zero, normalization skipped.\n";
  }

  for (int i = 0; i < n_cycles - 1; i++)
  {
    for (int k = 0; k < size; k++)
    {
      Ta_list.push_back(Ta_list.at(k));
      //	curr_time.push_back(curr_time.at(k));
    }
  }

  for (int k = 0; k < size * n_cycles; k++)
  {
    cout << k << " " << Ta_list.at(k) << endl;
  } */

  // exit(0);

  //  T = 0.9;
  dt = T / (size - 1);
  dt_mech = dt;
  cout << "size: " << size << endl;
  cout << "total time: " << T << endl;
  cout << "Mechanical time step: " << dt_mech << endl;
  // volume = arma::vec(Ta_list.size());
  // p_ven = arma::vec(Ta_list.size());
  // p_art = arma::vec(Ta_list.size());
  // p_LA = arma::vec(Ta_list.size());

  ephy.setup(monofile, cellmodel, odesolver, dt, T, pr, pa);
  ephy.update_matrix(false);
  ephy.init();
  // ephy.setup_types(typfile);

  elas.config(mshfile, parfile);
  elas.set_output_step(false);
  elas.init();
  elas.setup_data_writer(curr_time.size());
  // elas.set_output_step(false);

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

  /* for (int i =0; i<elas.get_mesh().get_n_elements();i++)
  {
    lat.push_back(0.136);
  }  */

  lat.set_size(elas.get_mesh().get_n_elements());
  has_eikonal = false; 
  pugi::xml_node element_data = doc.child("mesh").child("element_data");
  if(element_data)
  {
    for (pugi::xml_node elem = element_data.child("element"); elem; elem = elem.next_sibling("element"))
    {
        pugi::xml_node eikonal_p_elem = elem.child("eikonal");
        if(eikonal_p_elem)
        {
          has_eikonal = true; 
          int index; 
          index = elem.attribute("id").as_int();
          lat(index) = std::stod(eikonal_p_elem.text().as_string());
        }
        else
        {
          if(has_eikonal)
          {
            std::cerr << "ERROR: Some elements have a local activation time, and others do not.\n";
            exit(6);
          }
        }
    }
  }

  if(has_eikonal)
  {
    double min_val = lat.min(); 
    double max_val = lat.max(); 

    std::cout<<"Local activation time"<<std::endl;

    //TODO: read this information from file? 
    double begin_active_stress = 0.136;
    double latest_lat = 0.136 + 0.146;

    lat = begin_active_stress + (lat - min_val) * (latest_lat - begin_active_stress)/ (max_val-min_val);
    std::cout << " Earliest activation: " << lat.min() << "  Latest activation: " << lat.max()  << std::endl; 
    std::cout << " Loaded LATs: " << lat.n_elem << " values.\n";
  }
  else
  {
    lat.set_size(1);
    lat.fill(0.136); //if there isn't lat in the mesh file, we use 0.136 for all elements.
  }


  cout << "Precomputing active stress..." << endl;

  //TODO: define a generic active stress model 
  max_ta = 0; 
  for (int it_time = 0; it_time < curr_time.size(); it_time++)
  {
    double Ta_calculated =  solveTa((curr_time[it_time] / 1000.0) - 0.136, dt_mech); // using lat=0.136 only for instance, doesn't change the maximum value  
    if(Ta_calculated > max_ta)
        max_ta = Ta_calculated; 

  }

  cout << " Maximum active stress value (before rescaling): " << max_ta << std::endl;
  cout << " Maximum active stress value (after rescaling): "  << T_ref << std::endl;
  //obs: normalize the active stress with max_ta and multiply to T_ref. T_ref is aproximaly the maximum active stress value

  ta.zeros(nelem);
  dta.zeros(nelem);
  lc = 1.9;
  //TODO: Multiplos ciclos
}


void CardiacElectromechanic::Solve_System(double tt, double pressure, double pressure2)
{
  int size = elas.get_mesh().get_n_elements();
  int index = static_cast<int>(tt*1000); 

  if (has_eikonal)
  {
    for (int iel = 0; iel < size ; iel++)
    {
      double ta_calculated = solveTa(tt - lat(iel), dt_mech);
      ta_calculated *= (T_ref/max_ta);
      dta(iel) = ta_calculated - ta(iel);
      ta(iel) = ta_calculated; 
    }
  }
  else
  {
    double ta_calculated = solveTa(tt - lat(0), dt_mech);
    ta_calculated *= (T_ref/max_ta);
    dta.fill(ta_calculated - ta(0));
    ta.fill(ta_calculated);
  }

  P0 = pressure;
  
  cout << "Pressure: " << pressure << "Ta (mean): " << arma::mean(ta) << " dTA (mean): " <<arma::mean(dta) << " tt: " << tt << endl;
  cout << "Ta: min=" << ta.min() << " max=" << ta.max() 
      << " | dTa: min=" << dta.min() << " max=" << dta.max()
      << " | lat: min=" << lat.min() << " max=" << lat.max()
      << endl;
  elas.set_pressure_Ta(30, pressure, 20, pressure * 0.2, ta, dta);
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

  p_lv.push_back(0.0);
  p_rv.push_back(0.0);

  cout << "Solve 0 " << endl;
  Solve_System(0, p_lv[0], p_rv[0]);
  volume.push_back(elas.total_volume_cavity());

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
      cout << "Incrementa tempo..." << endl;
      tip.increase_time();
      cout << "\nTime: " << tip.time() << endl;

      i += 1;
      itempo = itempo + 1;

      double err = 1.;
      cout << "Acessa pressao anterior" << endl;
      p_0 = p_lv[i - 1];

      if (i >= 5)
      {
        cout << "Calcula Adams" << endl;
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

      Solve_System(tip.time(), p_0, p_0); 
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
        Solve_System(tip.time(), p_1, p_1);
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

      volume.push_back(elas.total_volume_cavity());

      timePoints.push_back(tip.time() + total_time);
      activeStressCurve.push_back(Ta0);

      cout << "Volume: " << volume[i] << endl;
      cout << "Pressure: " << p_lv[i] << endl;
      cout << "Ta: " << 0.0 << endl; //TODO: Ta value

      pv_file << tip.time() + total_time << " " << p_lv[i] << " " << volume[i] << " "
              << Ta0 << " "
              << p_art[i] << " " << p_ven[i] << " " << p_LA[i] << " "
              << V_art0 << " " << V_ven0 << " " << V_LA0 << " " << qao << " " << qmv << " " << qven << " " << qper << endl;

      cout << itempo << " " << tip.time() + total_time << " " << curr_time.at(itempo) << " " << p_lv[i] << " " << volume[i] << " "
           << Ta0 << " " << p_art[i] << " " << p_ven[i] << " " << p_LA[i] << " " << V_art0 << " " << V_ven0 << " " << V_LA0 << endl;

      if (k == 0)
      {
        ii++;
        cout << "Salvando XDMF... " << "Passo: " << ii << endl;
        elas.output_vtk(0, ii);
        elas.storeStress(ii);
      }
    }
  }

  pv_file.close();

  string stress_filename = filename + "_active_stress.txt";
  saveActiveStressToFile(stress_filename);

  elas.timer.summary();
  timer.summary();
}
