#include "cell_coupling.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "util/command_line_args.h"
#include "odes/torord_land.hpp"   // land_active_stiffness()

using std::cout;
using std::endl;

CellProfile CellProfile::for_model(const std::string & name, double T_ref)
{
  CellProfile p;

  if (name == "ToRORdLand")
  {
    // Ionico: Ta e uma variavel de estado (49), nao um valor monitorado, e a
    // celula precisa de fato despolarizar -> estimulo em cada no ao atingir o
    // LAT. Land entrega Ta em kPa; qualquer ganho adicional vem do ta_scale
    // por material do XML, para que um unico lugar o controle.
    p.ta_var  = 49;
    p.ta_unit = 1000.0;            // Land entrega kPa; a mecanica usa Pa
    p.vm_var  = 0;                 // Vm e a variavel de estado 0

    // Estados de distorcao de Land, so para o diagnostico de Ta negativa.
    // Indices conforme as leituras no topo de torord_land.cpp.
    p.land_xs = 43; p.land_xw = 44; p.land_zetas = 47; p.land_zetaw = 48;

    p.stim_amplitude = CommandLineArgs::read("-stimamp", -53.0);
    p.stim_duration  = CommandLineArgs::read("-stimdur", 2.0);   // ms
    p.uses_stretch   = true;

    // Euler explicito no ToRORd: 0.02 ms e o uso documentado deste solver. O
    // electromech usa 0.001 ms, 20x mais caro por passo de acoplamento e mais
    // seguro. Se 0.02 divergir, a rodada para com NaN em Ta e manda reduzir -dt.
    p.default_dt_ms  = 0.02;

    cout << " Active tension: state variable 49 (kPa) x 1000 -> Pa" << endl;
    cout << " LAT stimulus: amplitude " << p.stim_amplitude << ", duration "
         << p.stim_duration << " ms" << endl;
  }
  else if (name == "Kerkoff2003")
  {
    // Fenomenologico (Kerckhoffs): Ta e o valor monitorado 0, o modelo e
    // disparado pelo proprio LAT -> sem corrente de estimulo.
    p.ta_var  = -1;
    p.ta_unit = T_ref;
    p.default_dt_ms = 1.0;         // Euler explicito, escrito em s: 1 ms e estavel
  }
  else
  {
    throw std::runtime_error(
        "Modelo celular '" + name + "' sem mapeamento de tensao ativa no "
        "cardiac_cycle. Modelos suportados: Kerkoff2003, ToRORdLand. Para "
        "adicionar um, inclua um ramo em CellProfile::for_model() "
        "(cell_coupling.cpp) com o indice do estado de Ta e a conversao de "
        "unidade; electromechanic.cpp ja lista os indices de varios modelos.");
  }

  return p;
}

void LandCoupling::active_stiffness(const Cells & cells, const arma::vec & lambda,
                                    arma::vec & ka)
{
  const arma::uword n = lambda.n_elem;

  if (buf_xs_.n_elem != n)
  {
    buf_xs_.zeros(n);
    buf_xw_.zeros(n);
    buf_zetas_.zeros(n);
    buf_zetaw_.zeros(n);
  }

  cells.get_var(xs_, buf_xs_);
  cells.get_var(xw_, buf_xw_);

  ka.set_size(n);
  for (arma::uword i = 0; i < n; i++)
    ka(i) = land_active_stiffness(lambda(i), buf_xs_(i), buf_xw_(i));
}

void LandCoupling::report_negative_ta(const Cells & cells, const arma::vec & ta,
                                     const arma::vec & lambda_rate, double rate_clip)
{
  const int n_neg = static_cast<int>(arma::accu(ta < 0.0));

  cout << "  [MEF  ] *** WARNING: negative active tension at " << n_neg
       << "/" << ta.n_elem << " node(s), min = " << ta.min() << endl;

  // Without the Land state indices there is nothing further to say; the
  // decomposition below is specific to that submodel.
  if (xs_ < 0 || xw_ < 0
      || zetas_ < 0 || zetaw_ < 0)
    return;

  // Cells::get_var writes into v(i) without resizing it, exactly like the
  // pre-sized `ta` and `vm_node` buffers filled elsewhere. Sizing these here
  // rather than in config() keeps the allocation next to its only use and
  // costs nothing after the first call.
  if (buf_xs_.n_elem != ta.n_elem)
  {
    buf_xs_.zeros(ta.n_elem);
    buf_xw_.zeros(ta.n_elem);
    buf_zetas_.zeros(ta.n_elem);
    buf_zetaw_.zeros(ta.n_elem);
  }

  cells.get_var(xs_,    buf_xs_);
  cells.get_var(xw_,    buf_xw_);
  cells.get_var(zetas_, buf_zetas_);
  cells.get_var(zetaw_, buf_zetaw_);

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

  const bool have_rate = (lambda_rate.n_elem == ta.n_elem);

  for (arma::uword i = 0; i < ta.n_elem; i++)
  {
    if (ta(i) >= 0.0) continue;

    if (buf_zetas_(i) + 1.0 < 0.0) n_zetas++;
    else                           n_zetaw++;

    if (ta(i) < worst_ta)
    {
      worst_ta    = ta(i);
      xs_at_worst = buf_xs_(i);
      xw_at_worst = buf_xw_(i);
      worst_rate  = have_rate ? lambda_rate(i) : 0.0;
    }
  }

  cout << "  [MEF  ]   dominant term: ZETAW*XW at " << n_zetaw
       << " node(s), (ZETAS+1)<0 at " << n_zetas << " node(s)" << endl;
  cout << "  [MEF  ]   at the worst node: XS = " << xs_at_worst
       << ", XW = " << xw_at_worst
       << ", ZETAS = " << buf_zetas_.min()
       << ", ZETAW = " << buf_zetaw_.min() << endl;

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

    if (n_zetaw > n_zetas && rate_clip > 0.0
        && std::fabs(rate_crit) < rate_clip)
      cout << "  [MEF  ]   the threshold is BELOW -lamratemax ("
           << rate_clip << "): tightening the clip alone"
           << " will not fix this." << endl;
  }
}
// ======================================================================
//  Tensao ativa prescrita (-ta_curve)
// ======================================================================

void TaCurve::load(const std::string & path, double period_ms)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("-ta_curve: nao foi possivel abrir '" + path + "'");

  t_.clear();
  ta_.clear();
  std::string line;
  while (std::getline(in, line))
  {
    const std::size_t c = line.find('#');
    if (c != std::string::npos) line.erase(c);
    std::istringstream ss(line);
    double t, ta;
    if (!(ss >> t >> ta)) continue;          // linha vazia ou so comentario
    if (!t_.empty() && t <= t_.back())
      throw std::runtime_error("-ta_curve: o tempo tem de ser crescente em '" + path + "'");
    t_.push_back(t);
    ta_.push_back(ta);
  }
  if (t_.size() < 2)
    throw std::runtime_error("-ta_curve: '" + path + "' precisa de pelo menos 2 linhas (t [ms], Ta [Pa])");

  period_ = period_ms;

  cout << " Active tension: PRESCRIBED from '" << path << "' (" << t_.size()
       << " points, t = " << t_.front() << ".." << t_.back() << " ms, Ta max = "
       << *std::max_element(ta_.begin(), ta_.end()) << " Pa)" << endl;
  if (period_ > 0.0)
    cout << "   repeated every " << period_ << " ms; the cell model is still integrated"
         << " but its Ta is ignored" << endl;
}

double TaCurve::operator()(double t_ms) const
{
  double t = t_ms;
  if (period_ > 0.0)
    t = std::fmod(std::fmod(t, period_) + period_, period_);

  if (t <= t_.front()) return ta_.front();
  if (t >= t_.back())  return ta_.back();

  const std::size_t i = std::upper_bound(t_.begin(), t_.end(), t) - t_.begin();
  const double w = (t - t_[i - 1]) / (t_[i] - t_[i - 1]);
  return (1.0 - w) * ta_[i - 1] + w * ta_[i];
}
