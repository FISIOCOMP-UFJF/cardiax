#include "cellmodel.hpp"
#include "ode_solver.hpp"
#include "fitz_hugh_nagumo.hpp"
#include "nash_panfilov.hpp"
#include "mitchell_schaeffer.hpp"
#include "ten_tusscher2004.hpp"
#include "ten_tusscher2006.hpp"
#include "ten_tusscher_ta.hpp"
#include "rice_ord.hpp"
#include "rice_ten_tusscher.hpp"
#include "luo_rudy.hpp"
#include "mv.hpp"
#include "ord.hpp"
#include "simple_ode.hpp"
#include "ord.hpp"
#include "rice_ord.hpp"
#include "mynash.hpp"
#include "minimalmodel.hpp"
#include "mm_silva.hpp"
#include "torord_land.hpp"
#include "kerkoff2003.hpp"
#include "hdf5.h"

CellModel::CellModel(int num_states)
    : num_state_vars(num_states), f1(num_states), f2(num_states), monitored()
{
  // do nothing for now...
}

CellModel::~CellModel()
{
  delete ts;
  delete ode_solver;
}

void CellModel::advance(double *statevars, double t, double dt)
{
  i_stim = 0.0;
  ode_solver->advance(statevars, t, dt);
}

void CellModel::advance(double *statevars, double t, double dt, double istim)
{
  i_stim = istim;
  ode_solver->advance(statevars, t, dt);
}

void CellModel::compute_jacobian(double *states, double t, arma::mat &jac)
{
  int i, j;
  double max, yold, delta, temp;

  // Evaluate ODE rhs
  equation(t, states, f1.memptr());

  for (j = 0; j < num_state_vars; j++)
  {
    // Safe state variable i
    yold = states[j];

    // Step size
    max = 1e-5 > std::fabs(yold) ? 1e-5 : std::fabs(yold);
    delta = std::sqrt(1e-15 * max);

    // Trick
    temp = yold - delta;
    delta = temp - yold;

    // Pertubate state variables that is F(x + h ei)
    states[j] += delta;

    // Evaluate ODE rhs at F(x + h ei)
    equation(t, states, f2.memptr());

    // Approximates the j-th column
    for (i = 0; i < num_state_vars; i++)
      jac(i, j) = (f2(i) - f1(i)) / delta;

    // Restores the i-th value of the state variables
    states[j] = yold;
  }
}

void CellModel::set_celltype(int t)
{
  switch (t)
  {
  case 0:
    type = EPI;
    break;
  case 1:
    type = MCELL;
    break;
  case 2:
    type = ENDO;
    break;
  case 3:
    type = APEX;
    break;
  case 4:
    type = BASE;
    break;
  }
}

void CellModel::setup(string method, double timestep, double tf, double tp)
{
  ode_solver = ODESolver::create(method, this);
  ts = new TimeStepper(timestep, tf, tp);
}

void CellModel::solve()
{
  double t = 0.0;
  double dt = ts->timestep();
  ofstream out("ap");

  // Vector to hold variables
  arma::vec y(num_state_vars);

  // Set initial conditions
  init(y.memptr());

  // Write to file
  out << t << "\t" << y(0) << endl;

  // Start solving
  while (!(ts->finished()))
  {
    ts->increase_time();
    t = ts->time();

    // Write to file for test
    if (ts->time_to_print())
      out << (double)t << "\t" << y(0) << endl;

    // Advance in time using ODESolver
    ode_solver->advance(y.memptr(), t, dt);
  }

  out.close();
}

void CellModel::solveTest(double stim, double sstart, double sstop,
                          const std::string &fname)
{
  bool apply_stimulus;
  double t = 0.0;
  double dt = ts->timestep();
  ofstream out(fname.c_str());

  // Vector to hold variables
  arma::vec y(num_state_vars);

  // Set initial conditions
  init(y.memptr());

  out << fixed << setprecision(8) << t << "\t";
  for (int i = 0; i < num_state_vars; i++)
    out << scientific << setprecision(8) << y(i) << "\t";
  if (get_num_monitored() > 0)
  {
    for (int i = 0; i < get_num_monitored(); i++)
      out << scientific << setprecision(8) << get_monitored_value(i) << "\t";
  }
  out << endl;

  // Start solving
  while (!(ts->finished()))
  {
    ts->increase_time();
    t = ts->time();

    // Check and then apply stimulus if necessary
    apply_stimulus = (t >= sstart) && (t <= sstop);

    if (apply_stimulus)
      set_stimulus(stim);
    else
      set_stimulus(0.0);

    // Advance in time using ODESolver
    ode_solver->advance(y.memptr(), t, dt);

    if (ts->time_to_print())
    {
      // Write data
      out << fixed << setprecision(8) << t << "\t";

      for (int i = 0; i < num_state_vars; i++)
        out << scientific << setprecision(8) << y(i) << "\t";

      if (get_num_monitored() > 0)
      {
        for (int i = 0; i < get_num_monitored(); i++)
          out << scientific << setprecision(8) << get_monitored_value(i) << "\t";
      }
      out << endl;
    }
    // cout << " time " << t << " ms" << endl;
  }

  out.close();

  cout << " Done" << endl
       << endl;
}
/*
void CellModel::solveTestHDF5(double stim, double stime, double sdur,
                              const std::string &fname)
{
  bool apply_stimulus;
  double t = 0.0;
  double dt = ts->timestep();
  double sstart = stime;
  double sstop = stime + sdur;

  // prepare to write HDF5 file
  hsize_t nt = ((ts->endtime() - ts->starttime()) / ts->time_to_save()) + 1;
  hsize_t dims[1];
  hid_t file_id, group_id, dataset_id, dataspace_id;
  herr_t status;

  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  group_id = H5Gcreate2(file_id, "/protocol_1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // matrix to hold data
  int tam = ((ts->endtime() - ts->starttime()) / ts->time_to_save()) + 1;
  int size = num_state_vars + get_num_monitored() + 1;
  double **data_to_write = new double *[size];
  for (int i = 0; i < size; i++)
    data_to_write[i] = new double[tam];

  //
  // start simulating
  //

  // Vector to hold variables
  arma::vec y(num_state_vars);

  // Set initial conditions
  init(y.memptr());

  data_to_write[0][0] = t;
  int j = 0;
  for (int i = 0; i < size - get_num_monitored() - 1; i++)
    data_to_write[i + 1][0] = y(i);
  j++;

  if (get_num_monitored() > 0)
  {
    for (int i = j + 1; i < size; i++)
      data_to_write[i][0] = get_monitored_value(i);
  }

  // Start solving
  j = 0;
  int k = 1;
  while (!(ts->finished()))
  {
    ts->increase_time();
    t = ts->time();

    // Check and then apply stimulus if necessary
    apply_stimulus = (t >= sstart) && (t <= sstop);

    if (apply_stimulus)
      set_stimulus(stim);
    else
      set_stimulus(0.0);

    // Advance in time using ODESolver
    ode_solver->advance(y.memptr(), t, dt);

    if (ts->time_to_print())
    {
      data_to_write[0][k] = t;

      for (int i = 0; i < num_state_vars; i++)
        data_to_write[i + 1][k] = y(i);
      j++;

      if (get_num_monitored() > 0)
      {
        for (int i = j + 1; i < size; i++)
          data_to_write[i][k] = get_monitored_value(i);
      }

      k++;
    }
  }

  // write matrix's data to HDF5
  std::string dataset_name;
  std::map<int, std::string>::iterator it;
  it = var_names.begin();
  for (int i = 0; i < size; i++)
  {
    if (i == 0)
      dataset_name = "/protocol_1/t";
    else
    {
      dataset_name = "/protocol_1/";
      std::stringstream sstm;
      sstm << dataset_name << it->second;
      string aux = sstm.str();
      dataset_name = aux;
      ++it;
    }
    cout << "Writing dataset " << dataset_name << endl;

    // open HDF5 dataset and write data
    dims[0] = nt;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                           dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data_to_write[i]);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
  }

  status = H5Gclose(group_id);
  status = H5Fclose(file_id);

  if (status != 0)
    cout << "HDF5 Error at solveTestBCL" << endl;

  // free memory
  for (int i = 0; i < size; i++)
    delete[] data_to_write[i];
  delete[] data_to_write;

  cout << "Done" << endl;
}

void CellModel::solveTestBCL(double stim, double sdur, double bcl,
                             const std::string &fname)
{
  bool apply_stimulus;
  double t = 0.0;
  double dt = ts->timestep();
  double sstart = 1.0;
  double sstop = sstart + sdur;

  // prepare to write HDF5 file
  hsize_t nt = ((ts->endtime() - ts->starttime()) / ts->time_to_save()) + 1;
  hsize_t dims[1];
  hid_t file_id, group_id, dataset_id, dataspace_id;
  herr_t status;

  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  group_id = H5Gcreate2(file_id, "/protocol_1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // matrix to hold data
  int tam = ((ts->endtime() - ts->starttime()) / ts->time_to_save()) + 1;
  int size = num_state_vars + get_num_monitored() + 1;
  double **data_to_write = new double *[size];
  for (int i = 0; i < size; i++)
    data_to_write[i] = new double[tam];

  //
  // start simulating
  //

  // Vector to hold variables
  arma::vec y(num_state_vars);

  // Set initial conditions
  init(y.memptr());

  data_to_write[0][0] = t;
  int j = 0;
  for (int i = 0; i < size - get_num_monitored() - 1; i++)
    data_to_write[i + 1][0] = y(i);
  j++;

  if (get_num_monitored() > 0)
  {
    for (int i = j + 1; i < size; i++)
      data_to_write[i][0] = get_monitored_value(i);
  }

  // Start solving
  j = 0;
  int k = 1;
  while (!(ts->finished()))
  {
    ts->increase_time();
    t = ts->time();

    // Check and then apply stimulus if necessary
    apply_stimulus = (t >= sstart) && (t <= sstop);

    if (apply_stimulus)
    {
      set_stimulus(stim);

      if (t >= sstop)
      {
        sstart += bcl;
        sstop = sstart + sdur;
      }
    }
    else
      set_stimulus(0.0);

    // Advance in time using ODESolver
    ode_solver->advance(y.memptr(), t, dt);

    if (ts->time_to_print())
    {
      data_to_write[0][k] = t;

      for (int i = 0; i < num_state_vars; i++)
        data_to_write[i + 1][k] = y(i);
      j++;

      if (get_num_monitored() > 0)
      {
        for (int i = j + 1; i < size; i++)
          data_to_write[i][k] = get_monitored_value(i);
      }

      k++;
    }
  }

  // write matrix's data to HDF5
  std::string dataset_name;
  std::map<int, std::string>::iterator it;
  it = var_names.begin();
  for (int i = 0; i < size; i++)
  {
    if (i == 0)
      dataset_name = "/protocol_1/t";
    else
    {
      dataset_name = "/protocol_1/";
      std::stringstream sstm;
      sstm << dataset_name << it->second;
      string aux = sstm.str();
      dataset_name = aux;
      ++it;
    }
    cout << "Writing dataset " << dataset_name << endl;

    // open HDF5 dataset and write data
    dims[0] = nt;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                           dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data_to_write[i]);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
  }

  status = H5Gclose(group_id);
  status = H5Fclose(file_id);

  if (status != 0)
    cout << "HDF5 Error at solveTestHDF5" << endl;

  // free memory
  for (int i = 0; i < size; i++)
    delete[] data_to_write[i];
  delete[] data_to_write;

  cout << "Done" << endl;
}

void CellModel::solveTestRTT(double stim, double sdur, double itl, double bcl, double delta,
                             const std::string &fname)
{
  bool apply_stimulus;
  double t = 0.0;
  double dt = ts->timestep();
  double sstart = 1.0;
  double sstop = sstart + sdur;

  // prepare to write HDF5 file
  hsize_t nt = ((ts->endtime() - ts->starttime()) / ts->time_to_save()) + 1;
  hsize_t dims[1];
  hid_t file_id, group_id, dataset_id, dataspace_id;
  herr_t status;

  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  group_id = H5Gcreate2(file_id, "/protocol_1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // matrix to hold data
  int tam = ((ts->endtime() - ts->starttime()) / ts->time_to_save()) + 1;
  int size = num_state_vars + get_num_monitored() + 1;
  double **data_to_write = new double *[size];
  for (int i = 0; i < size; i++)
    data_to_write[i] = new double[tam];

  //
  // start simulating
  //

  // Vector to hold variables
  arma::vec y(num_state_vars);

  // Set initial conditions
  init(y.memptr());

  data_to_write[0][0] = t;
  int j = 0;
  for (int i = 0; i < size - get_num_monitored() - 1; i++)
    data_to_write[i + 1][0] = y(i);
  j++;

  if (get_num_monitored() > 0)
  {
    for (int i = j + 1; i < size; i++)
      data_to_write[i][0] = get_monitored_value(i);
  }

  // Start solving
  j = 0;
  int k = 1;
  while (!(ts->finished()))
  {
    ts->increase_time();
    t = ts->time();

    // Check and then apply stimulus if necessary
    apply_stimulus = (t >= sstart) && (t <= sstop);

    if (apply_stimulus)
    {
      set_stimulus(stim);

      if (t >= sstop && t < itl)
      {
        sstart += bcl;
        sstop = sstart + sdur;
      }
      if (t >= sstop && t >= itl)
      {
        bcl -= delta;
        sstart += bcl;
        sstop = sstart + sdur;
      }
    }
    else
      set_stimulus(0.0);

    // Advance in time using ODESolver
    ode_solver->advance(y.memptr(), t, dt);

    if (ts->time_to_print())
    {
      data_to_write[0][k] = t;

      for (int i = 0; i < num_state_vars; i++)
        data_to_write[i + 1][k] = y(i);
      j++;

      if (get_num_monitored() > 0)
      {
        for (int i = j + 1; i < size; i++)
          data_to_write[i][k] = get_monitored_value(i);
      }

      k++;
    }
  }

  // write matrix's data to HDF5
  std::string dataset_name;
  std::map<int, std::string>::iterator it;
  it = var_names.begin();
  for (int i = 0; i < size; i++)
  {
    if (i == 0)
      dataset_name = "/protocol_1/t";
    else
    {
      dataset_name = "/protocol_1/";
      std::stringstream sstm;
      sstm << dataset_name << it->second;
      string aux = sstm.str();
      dataset_name = aux;
      ++it;
    }
    cout << "Writing dataset " << dataset_name << endl;

    // open HDF5 dataset and write data
    dims[0] = nt;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                           dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data_to_write[i]);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
  }

  status = H5Gclose(group_id);
  status = H5Fclose(file_id);

  if (status != 0)
    cout << "HDF5 Error at solveTestRTT" << endl;

  // free memory
  for (int i = 0; i < size; i++)
    delete[] data_to_write[i];
  delete[] data_to_write;

  cout << "Done" << endl;
}
*/
CellModel *CellModel::create(std::string cellname)
{
  cout << "Cell model: " << cellname << endl;

  CellModel *ptr = NULL;

  if (cellname == "FHN")
    ptr = new FitzHughNagumo();
  else if (cellname == "NP")
    ptr = new NashPanfilov();
  else if (cellname == "MNP")
    ptr = new MyNashPanfilov();
  else if (cellname == "MINI")
    ptr = new MinimalModel();
  else if (cellname == "MMSilva")
    ptr = new MMSilva();
  else if (cellname == "MS")
    ptr = new MitchellSchaeffer();
  else if (cellname == "MV")
    ptr = new MinimalVentricular();
  else if (cellname == "LR1")
    ptr = new LuoRudy();
  else if (cellname == "ORd")
    ptr = new OHaraRudy();
  else if (cellname == "TNNP")
    ptr = new TenTusscher2004();
  else if (cellname == "TT2")
    ptr = new TenTusscher2006();
  else if (cellname == "TT2Ta")
    ptr = new TenTusscherTa();
  else if (cellname == "RiceTT2")
    ptr = new RiceTenTusscher();
  else if (cellname == "RiceORd")
    ptr = new RiceOHaraRudy();
  else if (cellname == "SODE")
    ptr = new SimpleODE();
  else if (cellname == "ToRORdLand")
    ptr = new TorordLand();
  else if(cellname == "Kerkoff2003")
    ptr = new Kerkoff2003(); 
  else
  {
    help(); 
    throw std::invalid_argument("Unknown CellModel.");
  }
  return ptr;
}

void CellModel::help()
{
    std::cout << "Available Cell Models:\n\n";

    std::cout << "  FHN          - FitzHugh-Nagumo\n";
    std::cout << "  NP           - Nash-Panfilov\n";
    std::cout << "  MNP          - My Nash-Panfilov\n";
    std::cout << "  MINI         - Minimal Model\n";
    std::cout << "  MMSilva      - MM Silva\n";
    std::cout << "  MS           - Mitchell-Schaeffer\n";
    std::cout << "  MV           - Minimal Ventricular\n";
    std::cout << "  LR1          - Luo-Rudy I\n";
    std::cout << "  ORd          - O'Hara-Rudy\n";
    std::cout << "  TNNP         - Ten Tusscher 2004\n";
    std::cout << "  TT2          - Ten Tusscher 2006\n";
    std::cout << "  TT2Ta        - Ten Tusscher (Ta)\n";
    std::cout << "  RiceTT2      - Rice + Ten Tusscher\n";
    std::cout << "  RiceORd      - Rice + O'Hara-Rudy\n";
    std::cout << "  SODE         - Simple ODE\n";
    std::cout << "  ToRORdLand   - ToR-ORd Land\n";
    std::cout << "  Kerkoff2003  - Kerkoff 2003\n";

    std::cout << "\nUsage example:\n";
    std::cout << "  -c ORd\n";
}
