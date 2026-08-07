#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <fstream>
#include "cells.hpp"
#include "cellmodel.hpp"
#include "ode_solver.hpp"
#include "fitz_hugh_nagumo.hpp"
#include "ten_tusscher2006.hpp"
#include "ten_tusscher_ta.hpp"
#include "rice_ten_tusscher.hpp"
#include "luo_rudy.hpp"
#include "mv.hpp"
#include "torord_land.hpp"
#include "util/command_line_args.h"


using namespace std;

// FHN 0.1   20.0
// LRI 0.01  500.0
// TT2 0.001 500.0

void usage()
{
  cout << "Usage: \n";
  cout << " ./testCellmodel -c <Model> -m <Method> -dt <Step> -p <Protocol> [More options]\n" << endl;
  cout << "Options:" << endl;
  cout << " -p protocol (single, pacing, restitution)" << endl;
  cout << " -c model" << endl;
  cout << " -m method (ExplicitEuler, ImplicitEuler, RungeKutta4)" << endl;
  cout << " -T total time" << endl;
  cout << " -dt step" << endl;
  cout << " -st stimulus time" << endl;
  cout << " -sv stimulus value" << endl;
  cout << " -sd stimulus duration" << endl;
  cout << " -tr time rate" << endl;
  cout << " -ct cell type (EPI,MCELL,ENDO,APEX,BASE - TT2 and TT2Ta onyl)" << endl;
  cout << " -out output file name" << endl;
  cout << " -bcl basic cycle length" << endl;
  cout << " -ncl number of cycles" << endl;
  cout << " -nrc number of restitution cycles" << endl;
  cout << " -d delta - time step (variation) for restitution protocol" << endl << endl;
  cout << " Models" << endl;
  cout << "   NP     : Nash-Panfilov" << endl;
  cout << "   MNP    : MyNash-Panfilov" << endl;
  cout << "   MS     : Mitchell-Schaeffer" << endl;
  cout << "   MV     : Minimal Ventricular" << endl;
  cout << "   FHN    : FitzHugh-Nagumo" << endl;
  cout << "   LR1    : Luo-Rudy I" << endl;
  cout << "   TT2    : ten Tusscher" << endl;
  cout << "   TT2Ta  : ten Tusscher + Active tension (Ta)" << endl;
  cout << "   RiceTT2: Rice + ten Tusscher" << endl;
  cout << "   SODE   : simple ODE for test" << endl;
  cout << endl;
}

int main(int argc, const char **argv)
{
//  int celltp = 0;
//  double step, T, stim, sini, send, tts;
//  string model;
//  string method;
//  string outf;
//
//
//  // TODO: include BCL protocol
//  // BCL1, BCL2, etc....
//
//  if (argc < 9)
//  {
//    cout << "\n Usage: ";
//    cout << "testCellmodel <Model> <Method> <Step> <T> <Sval> <Sini> <Send> <tts> [Type] [Output]\n" << endl;
//    cout << " Models" << endl;
//    cout << "   NP     : Nash-Panfilov" << endl;
//    cout << "   MNP    : MyNash-Panfilov" << endl;
//    cout << "   MS     : Mitchell-Schaeffer" << endl;
//    cout << "   MV     : Minimal Ventricular" << endl;
//    cout << "   FHN    : FitzHugh-Nagumo" << endl;
//    cout << "   LR1    : Luo-Rudy I" << endl;
//    cout << "   TT2    : ten Tusscher" << endl;
//    cout << "   TT2Ta  : ten Tusscher + Active tension (Ta)" << endl;
//    cout << "   RiceTT2: Rice + ten Tusscher" << endl;
//    cout << "   SODE   : simple ODE for test" << endl;
//    cout << endl;
//    cout << " Methods" << endl;
//    cout << "   ExplicitEuler" << endl;
//    cout << "   ImplicitEuler" << endl;
//    cout << "   RungeKutta4" << endl;
//    cout << endl;
//    cout << " Cell Type (TT2 and TT2Ta only)" << endl;
//    cout << "   0 : EPI" << endl;
//    cout << "   1 : MCELL" << endl;
//    cout << "   2 : ENDO" << endl;
//    cout << "   3 : APEX" << endl;
//    cout << "   4 : BASE" << endl;
//    cout << endl;
//    exit(1);
//  }
//
//  model  = argv[1];
//  method = argv[2];
//  step   = atof(argv[3]);
//  T      = atof(argv[4]);
//  stim   = atof(argv[5]);
//  sini   = atof(argv[6]);
//  send   = atof(argv[7]);
//  tts    = atof(argv[8]);
//
//  CellModel * cell = CellModel::create(model);
//  cell->setup(method, step, T, tts);
//
//  if(argc >= 10)
//  {
//    celltp = atoi(argv[9]);
//    cout << " Using celltype: " << celltp << endl;
//    if(model == "TT2Ta" || model == "RiceTT2" || model == "TNNP" ||
//       model == "ORd"   || model == "RiceORd" || model == "MV")
//      cell->set_celltype(celltp);
//  }
//
//  if(argc == 11)
//  {
//    outf = argv[10];
//    cell->solveTest(stim, sini, send, outf);
//  }
//  else
//  {
//    //cell->solveTest(stim, sini, send, "output.txt");
//    cell->solveTestHDF5(stim, sini, send, "output.h5");
//  }
//
//
//  delete cell;

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
