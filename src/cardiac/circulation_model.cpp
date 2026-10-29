#include "circulation_model.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int CirculationModel::limit_beats(const TimeBase & tb, int beats,
                                  double period_ms, double lead_in_ms)
{
  if (beats < 1) beats = 1;
  if (tb.ceiling_ms <= 0.0) return beats;

  const int fits = static_cast<int>(
      std::floor((tb.ceiling_ms - lead_in_ms) / period_ms + 1e-9));

  if (fits < beats)
  {
    std::cout << " *** -t " << tb.ceiling_ms << " ms cobre apenas " << fits
              << " batimento(s) de " << period_ms << " ms";
    if (lead_in_ms > 0.0)
      std::cout << " apos a fase 0 de " << lead_in_ms << " ms";
    std::cout << "; -beats " << beats << " reduzido para "
              << std::max(1, fits) << std::endl;
    beats = std::max(1, fits);
  }
  return beats;
}

void CirculationModel::print_status(const StepContext & s,
                                    const std::string & tail)
{
  // std::defaultfloat restaura o FORMATO mas nao a PRECISAO: sem salvar e
  // restaurar, o setprecision(1) vazaria para todas as linhas seguintes.
  const std::streamsize prec0 = std::cout.precision();
  std::cout << "\n[t = " << std::fixed << std::setprecision(1) << s.t_ms
            << " ms | step " << s.index + 1 << "/" << s.nsteps;
  if (!tail.empty()) std::cout << " | " << tail;
  std::cout << "]" << std::defaultfloat << std::endl;
  std::cout.precision(prec0);
}

void CirculationModel::open_history(const std::string & path,
                                    const std::string & header)
{
  history_.open(path.c_str());
  if (history_.is_open() && !header.empty())
    history_ << header << "\n";
}
const std::vector<double> & CirculationModel::state_entry(const State & st,
                                                          const std::string & key)
{
  State::const_iterator it = st.find(key);
  if (it == st.end())
    throw std::runtime_error("Checkpoint: estado da circulacao sem a entrada '"
                             + key + "'");
  return it->second;
}
