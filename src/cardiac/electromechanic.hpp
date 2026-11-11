#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include <cmath>
#include <fstream>
#include <string>
#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include "eikonal_reactive.hpp"

class Electromechanic
{
public:
  explicit Electromechanic(const std::string & epmodel) : ep_type(epmodel) {}
  ~Electromechanic() { delete ephy; }

  void config(const std::string & basename);
  void solve();
  CardiacProblem & ref() { return *ephy; }

  //! Checkpoint interval in ms; call after config() (uses dt_ephy)
  void set_checkpoint_interval(double interval_ms)
  {
    checkpoint_interval = (interval_ms > 0.0)
                        ? static_cast<int>(std::lround(interval_ms / dt_ephy)) : 0;
  }

  //! Print the command-line options of the coupled solver
  static void usage();

private:
  // --- setup, called from config() ---
  void read_options();              // command line + time-step checks
  void setup_electrophysiology();   // ephy: create, setup, init
  void setup_mechanics();           // elas: config, init, fibres
  void setup_output();              // fields, file, static data

  // --- models ---
  CardiacProblem * ephy = nullptr;  // Monodomain or MonodomainDeformation
  UpdatedLagrangian elas;

  // --- options ---
  std::string filename;
  std::string ep_type;
  std::string cell_model;
  std::string odesolver;

  // --- time ---
  double t_final = 0.0;
  double dt_ephy = 0.01;
  double dt_mech = 0.1;
  double dt_out  = 1.0;
  int steps_per_mech = 1;           // EP steps between mechanics solves
  int steps_per_out  = 1;           // EP steps between output frames
  int checkpoint_interval = 0;      // in EP steps; 0 = off

  // --- coupled fields ---
  arma::vec vm, ta, u_field;

  TimerSection timer;
};

#endif