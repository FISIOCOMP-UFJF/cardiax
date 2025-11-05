#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include <fstream>

class CardiacElectromechanic
{
public:
  CardiacElectromechanic(const std::string &epmodel);

  ~CardiacElectromechanic();

  void config(const string &basename);

  MonodomainDeformation &ref() { return ephy; }

  void solve();

public:
  double solveTa(double cur_time, double dt)
  {
    lc = lc + dt * (Ea * (ls0 - lc) - 1.0) * v0;
    return sigma(cur_time, ls0, lc);
  }

  double fiso(double lc)
  {
    return (lc <= a7) ? 0.0 : T0 * pow(tanh(a6 * (lc - a7)), 2);
  }

  double ftwich(double cur_time, double ls)
  {
    double tmax = b * (ls - ld);
    if (cur_time < 0.0 || cur_time > tmax)
      return 0.0;
      
    return pow(tanh(cur_time / tr), 2) * pow(tanh((tmax - cur_time) / td), 2);
  }

  double sigma(double cur_time, double ls, double lc)
  {
    return (ls / ls0) * fiso(lc) * ftwich(cur_time, ls) * (ls - lc) * Ea;
  }

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
  MonodomainDeformation ephy;
  UpdatedLagrangian elas;
  std::vector<arma::mat33 *> vec_stress;
  std::vector<arma::mat33 *> vec_fib;
  std::vector<arma::mat33 *> vec_fib0;
  TimerSection timer;
  std::vector<double> curr_time;
  std::vector<double> recorded_time;
  std::vector<double> p_lv, p_rv, Ta_list, volume, p_art, p_ven, p_LA;
  
  arma::mat Ta_matrix; 
  arma::vec dta;
  arma::vec ta;
  arma::vec lat;

  void Solve_System(double tt, double active_stress, double pressure, double pressure2);

  const double a6 = 2.0;
  const double a7 = 1.5;
  const double T0 = 180.0;
  const double Ea = 20.0;
  const double v0 = 7.5;
  const double ls0 = 1.9;
  const double tr = 0.075;
  const double td = 0.075;
  const double b = 0.21; 
  const double ld = -0.4;
  double lc = 0.0;
  double max_ta = 0.0; 

  std::vector<double> activeStressCurve;
  std::vector<double> timePoints;
  bool curveCalculated = false;
};

#endif
