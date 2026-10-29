#include "legacy_circulation.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "util/pugixml.hpp"

using std::cout;
using std::endl;

void LegacyCirculation::configure(const std::string & mshfile, const TimeBase & tb)
{
  base_ = mshfile;

  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());
  if (!result)
  {
    cout << "XML parsed with errors\n";
    cout << "Error description: " << result.description() << "\n";
    cout << "Error offset: " << result.offset << "\n\n";
  }

  pugi::xml_node pv = doc.child("mesh").child("pvloop");

  c_art_ = pv.attribute("C_art").as_double();
  c_ven_ = pv.attribute("C_ven").as_double();
  part_  = pv.attribute("p_art").as_double();
  pven_  = pv.attribute("p_ven").as_double();
  Rmv_   = pv.attribute("R_mv").as_double();
  Rven_  = pv.attribute("R_ven").as_double();
  Rper_  = pv.attribute("R_per").as_double();
  Rao_   = pv.attribute("R_ao").as_double();
  V_art_zero_ = pv.attribute("V_art_zero").as_double();
  V_ven_zero_ = pv.attribute("V_ven_zero").as_double();
  E_es_LA_ = pv.attribute("E_es_LA").as_double();
  A_LA_    = pv.attribute("A_LA").as_double();
  B_LA_    = pv.attribute("B_LA").as_double();
  Tmax_    = pv.attribute("T_max").as_double();
  tau_     = pv.attribute("tau").as_double();
  const bool biv = pv.attribute("biv").as_int() != 0;

  // Curva de pressao: um <step> por instante. Com biv, VE e VD tem curvas
  // proprias; senao so ha a do VE e o VD fica a zero.
  for (pugi::xml_node node = pv.child("step"); node; node = node.next_sibling("step"))
  {
    curr_time_.push_back(node.attribute("time").as_double());
    p_lv_.push_back(node.attribute(biv ? "pressure_lv" : "pressure").as_double());
    p_rv_.push_back(biv ? node.attribute("pressure_rv").as_double() : 0.0);
  }
  if (p_lv_.empty())
    throw std::runtime_error("Caminho legado: o <pvloop> da malha nao tem "
                             "nenhum <step>. Use -circ 1 ou -pw 1, ou forneca a curva.");
  cout << " pvloop curve: " << p_lv_.size() << " step(s)" << endl;

  // O <pvloop> usa as unidades do cardiax: total_time e <step time=...> em ms,
  // pressoes em Pa. O passo de acoplamento e a amostragem da curva: p_lv_[i] so
  // faz sentido se um passo avanca exatamente um <step>, entao -dt_mech nao vale
  // aqui. A EP pode usar um passo menor (-dt): o laco faz os sub-passos.
  const double total_ms = pv.attribute("total_time").as_double();
  const int    size     = static_cast<int>(p_lv_.size());
  if (size < 2 || total_ms <= 0.0)
    throw std::runtime_error("Caminho legado: <pvloop> sem total_time/steps validos.");
  if (total_ms < 50.0)
    throw std::runtime_error("<pvloop total_time=\"" + std::to_string(total_ms) +
        "\"> parece estar em SEGUNDOS. O cardiax usa ms: multiplique por 1000 "
        "total_time e todos os <step time=...>.");
  const double step_ms = total_ms / (size - 1);

  // Duracao: a da curva, ou a de -t (aqui -t e a duracao, nao um teto)
  double T_ms = total_ms;
  if (tb.ceiling_ms > 0.0 && tb.ceiling_ms != T_ms)
  {
    cout << " total_time from the mesh (" << T_ms << " ms) overridden by -t "
         << tb.ceiling_ms << " ms" << endl;
    T_ms = tb.ceiling_ms;
  }

  plan_.beats             = 1;
  plan_.run_ms            = T_ms;
  plan_.step_ms           = step_ms;
  plan_.nsteps            = static_cast<int>(std::round(T_ms / step_ms));
  plan_.ep_dt_hint_ms     = step_ms;      // sem -dt: um passo da EP por passo da curva
  plan_.period_ms         = 0.0;          // o warm up sai de T
  plan_.span_margin_steps = 0;
  plan_.default_out_every = 1;            // um frame por passo, como sempre
  plan_.store_stress      = true;
}

void LegacyCirculation::do_begin()
{
  // O vetor ja vem com a curva; um elemento extra guarda a ultima escrita.
  p_lv_.push_back(0.0);
  p_rv_.push_back(0.0);

  cout << "Solve 0 " << endl;
  mech_->apply_pressures(0.0, p_lv_[0], p_rv_[0]);
  volume_.push_back(mech_->volume_LV());
  volume_rv_.push_back(mech_->volume_RV());

  cout << "Volume: " << volume_[0] << endl;
  cout << "Pressure: " << p_lv_[0] << " " << p_rv_[0] << endl;

  p_ven_.push_back(pven_);
  p_art_.push_back(part_);
  p_LA_.push_back(0.0);

  V_LP0_  = volume_[0];
  V_art0_ = p_art_[0] * c_art_ + V_art_zero_;
  V_ven0_ = p_ven_[0] * c_ven_ + V_ven_zero_;
  V_LA0_  = 10.0;

  open_history(base_ + "_pvloop.txt", "");
}

void LegacyCirculation::announce(const StepContext & s)
{
  if (s.log) print_status(s, "");
}

void LegacyCirculation::advance(const StepContext & s)
{
  const double pi = 3.14159265358979323846;
  const double DT = s.dt_ms;
  const int    i  = s.index + 1;     // posicao deste passo nos historicos

  // ---- palpite inicial da pressao: Adams-Bashforth de 4a ordem ---------
  double p_0 = p_lv_[i - 1];
  double p_1;
  if (i >= 5)
  {
    double dpn  = p_lv_[i - 1] - p_lv_[i - 2];
    double dpn1 = p_lv_[i - 2] - p_lv_[i - 3];
    double dpn2 = p_lv_[i - 3] - p_lv_[i - 4];
    double dpn3 = p_lv_[i - 4] - p_lv_[i - 5];
    p_1 = p_0 + (DT / 24.) * (55. * dpn / DT - 59. * dpn1 / DT + 37. * dpn2 / DT - 9. * dpn3 / DT);
  }
  else
  {
    float DP = p_ven_[0];
    float DV = 251. - 199.;
    float Q = DV / 200.;
    p_1 = p_0 + (DP / DV) * Q * DT;
  }

  mech_->apply_pressures(s.t_ms, p_0, p_rv_[i - 1]);
  double Vf_0 = mech_->volume_total();
  double Vf_1 = Vf_0;

  double p_ven1 = p_ven_[i - 1];
  double p_art1 = p_art_[i - 1];
  double p_LA1  = p_LA_[i - 1];

  double V_LA = V_LA0_, V_LP = V_LP0_, V_art = V_art0_, V_ven = V_ven0_;
  double qmv = 0, qao = 0, qven = 0, qper = 0;

  // ---- Newton sobre o volume da cavidade ------------------------------
  double err = 1.;
  int iterations = 0;
  while (err > 0.001)
  {
    iterations += 1;

    mech_->apply_pressures(s.t_ms, p_1, p_rv_[i]);
    Vf_1 = mech_->volume_total();
    double C = (p_1 - p_0) / (Vf_1 - Vf_0);

    qmv  = (p_LA1 >= p_1)   ? (0.5 / Rmv_) * (p_LA1 + p_LA_[i - 1] - p_1 - p_lv_[i - 1]) : 0.0;
    qao  = (p_1 >= p_art1)  ? (0.5 / Rao_) * (p_1 + p_lv_[i - 1] - p_art1 - p_art_[i - 1]) : 0.0;
    qper = (0.5 / Rper_) * (p_art1 + p_art_[i - 1] - p_ven1 - p_ven_[i - 1]);
    qven = (0.5 / Rven_) * (p_ven1 + p_ven_[i - 1] - p_LA1 - p_LA_[i - 1]);

    V_LA  = V_LA0_  + (qven - qmv) * DT;
    V_LP  = V_LP0_  + (qmv - qao) * DT;
    V_art = V_art0_ + (qao - qper) * DT;
    V_ven = V_ven0_ + (qper - qven) * DT;

    p_art1 = (V_art - V_art_zero_) / c_art_;
    p_ven1 = (V_ven - V_ven_zero_) / c_ven_;

    // atrio esquerdo: elastancia variavel no tempo
    const double tc = curr_time_.at(i);
    double pES_LA = E_es_LA_ * (V_LA - V_LA_zero_);
    double pED_LA = A_LA_ * (exp(B_LA_ * (V_LA - V_LA_zero_)) - 1.0);
    double e_func = (tc >= 0 && tc <= (3. / 2.) * Tmax_)
                  ? 0.5 * (sin((pi / Tmax_) * tc - pi / 2.) + 1.0)
                  : 0.5 * exp(-(tc - (3. / 2.) * Tmax_) / tau_);
    p_LA1 = e_func * pES_LA + (1.0 - e_func) * pED_LA;

    double r_est = Vf_1 - V_LP;

    p_0  = p_1;
    Vf_0 = Vf_1;
    p_1  = p_1 - C * r_est;

    err = fabs((Vf_1 - V_LP) / V_LP);
    if (s.log)
      cout << "  [NEWTON] it " << iterations << ": error = " << err << endl;
  }

  p_art_.push_back(p_art1);
  p_ven_.push_back(p_ven1);
  p_LA_.push_back(p_LA1);
  V_LP0_  = V_LP;
  V_art0_ = V_art;
  V_ven0_ = V_ven;
  V_LA0_  = V_LA;

  p_lv_[i] = p_1;
  volume_.push_back(mech_->volume_LV());
  volume_rv_.push_back(mech_->volume_RV());

  if (s.log)
    cout << "  [STATE] pressure converged in " << iterations
         << " iteration(s): p_LV = " << p_lv_[i] << ", V_LV = " << volume_[i]
         << ", p_art = " << p_art_[i] << ", p_ven = " << p_ven_[i]
         << ", p_LA = " << p_LA_[i] << endl;

  if (history_.is_open())
    history_ << s.t_ms << " " << p_lv_[i] << " " << volume_[i] << " "
             << 0.0 << " "
             << p_art_[i] << " " << p_ven_[i] << " " << p_LA_[i] << " "
             << V_art0_ << " " << V_ven0_ << " " << V_LA0_ << " "
             << qao << " " << qmv << " " << qven << " " << qper << " "
             << p_rv_[i] << " " << volume_rv_[i] << endl;
}

void LegacyCirculation::end()
{
  if (history_.is_open()) history_.close();
}
void LegacyCirculation::save_state(State & st) const
{
  // Historicos inteiros: o Adams-Bashforth e as vazoes leem passos anteriores,
  // e p_lv_ guarda tambem o resto da curva tabelada.
  st["p_lv"]      = p_lv_;
  st["p_rv"]      = p_rv_;
  st["p_art"]     = p_art_;
  st["p_ven"]     = p_ven_;
  st["p_LA"]      = p_LA_;
  st["volume"]    = volume_;
  st["volume_rv"] = volume_rv_;
  st["reservoirs"] = { V_LP0_, V_art0_, V_ven0_, V_LA0_ };
}

void LegacyCirculation::restore_state(const State & st)
{
  p_lv_      = state_entry(st, "p_lv");
  p_rv_      = state_entry(st, "p_rv");
  p_art_     = state_entry(st, "p_art");
  p_ven_     = state_entry(st, "p_ven");
  p_LA_      = state_entry(st, "p_LA");
  volume_    = state_entry(st, "volume");
  volume_rv_ = state_entry(st, "volume_rv");

  const std::vector<double> & r = state_entry(st, "reservoirs");
  V_LP0_ = r.at(0); V_art0_ = r.at(1); V_ven0_ = r.at(2); V_LA0_ = r.at(3);
}
