#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include "eikonal.hpp"
#include <fstream>

class CardiacElectromechanic
{
public:
  CardiacElectromechanic(const std::string &epmodel);

  ~CardiacElectromechanic();

  void config(const string &basename);

  // MonodomainDeformation &ref() { return ephy; }
  Eikonal &ref() { return ephy; }
  void solve();

public:
  void saveActiveStressToFile(const std::string &filename)
  {
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
      std::cerr << "Erro: Nao foi possivel abrir o arquivo" << filename << std::endl;
    }
    for (size_t i = 0; i < activeStressCurve.size(); ++i)
    {
      outFile << timePoints[i] << " " << activeStressCurve[i] << "\n";
    }
    outFile.close();
    std::cout << "Dados da tensão ativa salvos em: " << filename << std::endl;
  }

private:
  double T_ref;
  double dt_mech;
  int n_cycles = 1;
  double Ta0 = 0.0;
  double P0;
  double dTa0;
  double c_art, c_ven, Rmv, Rven, Rao, Rper;
  double V_art_zero, V_ven_zero;
  double E_es_LA, A_LA, B_LA, Tmax, tau;
  double P_o, part, pven, stroke_volume;
  string filename;
  // MonodomainDeformation ephy;
  Eikonal ephy; 
  UpdatedLagrangian elas;
  std::vector<arma::mat33 *> vec_stress;
  std::vector<arma::mat33 *> vec_fib;
  std::vector<arma::mat33 *> vec_fib0;
  TimerSection timer;
  std::vector<double> curr_time;
  std::vector<double> recorded_time;
  std::vector<double> p_lv, p_rv, Ta_list, volume, p_art, p_ven, p_LA;
  
  arma::vec dta;
  arma::vec ta;
  arma::vec lat;

  void Solve_System(double tt, double pressure, double pressure2);


  std::vector<double> activeStressCurve;
  std::vector<double> timePoints;
  bool curveCalculated = false;
  bool has_eikonal = false; 
};

#endif
