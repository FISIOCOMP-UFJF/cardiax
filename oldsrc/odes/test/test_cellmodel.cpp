#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <fstream>
#include "../odes.h"
#include "util/command_line_args.h"

using namespace std;

void usage()
{
  const int W = 26; // column width for alignment
  auto opt = [&](const string& flag, const string& desc) {
    cout << "  " << left << setw(W) << flag << desc << "\n";
  };

  cout << "\n";
  cout << "  Cardiax - Cell Model Solver\n";
  cout << "  ============================\n\n";

  cout << "  Usage:\n";
  cout << "    ./testCellmodel -c <model> -m <method> -dt <step> -p <protocol> [options]\n\n";

  cout << "  Protocols (-p):\n";
  opt("-p single",       "Single stimulus run (default)");
  opt("-p pacing",       "Basic cycle length pacing");
  opt("-p restitution",  "Restitution protocol");
  cout << "\n";

  cout << "  Required options:\n";
  opt("-c  <model>",     "Cell model (see list below)");
  opt("-m  <method>",    "ODE method (see list below)");
  opt("-dt <step>",      "Time step (ms)");
  cout << "\n";

  cout << "  General options:\n";
  opt("-T   <time>",     "Total simulation time (ms)");
  opt("-tr  <rate>",     "Time rate");
  opt("-out <file>",     "Output file name [default: output.h5]");
  opt("-ct  <type>",     "Cell type: EPI | MCELL | ENDO | APEX | BASE");
  cout << "\n";

  cout << "  Stimulus options:\n";
  opt("-st  <time>",     "Stimulus start time (ms)");
  opt("-sv  <value>",    "Stimulus value");
  opt("-sd  <duration>", "Stimulus duration (ms) [default: 1.0]");
  cout << "\n";

  cout << "  Pacing options:\n";
  opt("-bcl <length>",   "Basic cycle length (ms)");
  opt("-nbc <n>",        "Number of basic cycles");
  cout << "\n";

  cout << "  Restitution options:\n";
  opt("-bcl <length>",   "Basic cycle length (ms)");
  opt("-nbc <n>",        "Number of basic cycles (to reach steady state)");
  opt("-nrc <n>",        "Number of restitution cycles");
  opt("-d   <delta>",    "Time step variation per restitution cycle (ms)");
  cout << "\n";

  cout << "  Models:\n";
  opt("NP",          "Nash-Panfilov");
  opt("MNP",         "MyNash-Panfilov");
  opt("MS",          "Mitchell-Schaeffer");
  opt("MV",          "Minimal Ventricular");
  opt("FHN",         "FitzHugh-Nagumo");
  opt("LR1",         "Luo-Rudy I");
  opt("TT2",         "ten Tusscher 2006");
  opt("TT2Ta",       "ten Tusscher + Active tension (cell types: EPI,MCELL,ENDO,APEX,BASE)");
  opt("RiceTT2",     "Rice + ten Tusscher");
  opt("ToRORdLand",  "ToR-ORd + Land");
  opt("SODE",        "Simple ODE (test)");
  cout << "\n";

  cout << "  Methods:\n";
  opt("ExplicitEuler",  "Forward Euler");
  opt("ImplicitEuler",  "Backward Euler");
  opt("RungeKutta4",    "4th order Runge-Kutta");
  cout << "\n";

  cout << "  Examples:\n";
  cout << "    ./testCellmodel -p single      -c TT2 -m RungeKutta4 -dt 0.001 -T 500 -sv -52 -st 10 -sd 1\n";
  cout << "    ./testCellmodel -p pacing      -c FHN -m ExplicitEuler -dt 0.1 -bcl 500 -nbc 10 -sv -0.3 -sd 1\n";
  cout << "    ./testCellmodel -p restitution -c TT2 -m RungeKutta4 -dt 0.001 -bcl 500 -nbc 20 -nrc 10 -d 10\n";
  cout << "\n";
}

int main(int argc, const char **argv)
{
  if(argc < 6)  // not enough params
  {
    usage();
    exit(1);
  }

  int nbc, nrc, celltp;
  double step, T, bcl, delta, stim, stime, sdur, tts;
  string protocol, model, method, outf;

  CommandLineArgs::init(argc, argv);

  protocol = CommandLineArgs::read("-p", "single");   // protocol
  model = CommandLineArgs::read("-c", model);         // model
  method = CommandLineArgs::read("-m", method);       // method
  step = CommandLineArgs::read("-dt", step);          // time step
  stim = CommandLineArgs::read("-sv", stim);          // stimulus value
  sdur = CommandLineArgs::read("-sd", 1.0);          // stimulus duration
  tts = CommandLineArgs::read("-tr", 1.0);            // time rate
  celltp = CommandLineArgs::read("-ct", 0);           // cell type
  outf = CommandLineArgs::read("-out", "output.h5");  // output file name

  CellModel * cell = CellModel::create(model);
  if(celltp != 0)
  {
    cout << " Using celltype: " << celltp << endl;
    if(model == "TT2Ta" || model == "RiceTT2" || model == "TNNP" ||
       model == "ORd"   || model == "RiceORd" || model == "MV" || model == "ToRORdLand")
      cell->set_celltype(celltp);
  }

  if(protocol == "single")
  {
    T = CommandLineArgs::read("-T", T);
    stime = CommandLineArgs::read("-st", stime);

    cell->setup(method, step, T, tts);
//    cell->solveTestHDF5(stim, stime, sdur, outf);
    cell->solveTest(stim, stime, sdur, outf);
  }
  else if(protocol == "pacing")
  {
    bcl = CommandLineArgs::read("-bcl", bcl);    // basic cycle length
    nbc = CommandLineArgs::read("-nbc", nbc);    // number of basic cycles
    T = bcl * nbc;

    cell->setup(method, step, T, tts);
    //cell->solveTestBCL(stim, sdur, bcl, outf);
  }
  else if(protocol == "restitution")
  {
    bcl = CommandLineArgs::read("-bcl", bcl);    // basic cycle length
    nbc = CommandLineArgs::read("-nbc", nbc);    // number of basic cycles
    nrc = CommandLineArgs::read("-nrc", nrc);    // number of restitution cycles
    delta = CommandLineArgs::read("-d", delta);  // time variation
    double itl = (bcl * nbc);                    // initial time loop, to reach stationary state
    T = itl;
    for(int i = 1; i <= nrc; i++)
      T = T + (bcl - (i*delta));

    cell->setup(method, step, T, tts);
    //cell->solveTestRTT(stim, sdur, itl, bcl, delta, outf);
  }

  delete cell;

}
