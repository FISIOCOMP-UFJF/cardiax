#include "conductivity.hpp"

#include <unordered_map>
#include <stdexcept>

const ConductivityValues & ConductivityModel::values(int region_id) const
{
  auto it = table_.find(region_id);
  return (it != table_.end()) ? it->second : default_;
}

arma::mat33 ConductivityModel::at(int region_id, int ndim,
                                  const arma::vec3 & f,
                                  const arma::vec3 & s,
                                  const arma::vec3 & n) const
{
  return assemble(ndim, values(region_id), f, s, n);
}

arma::mat33 ConductivityModel::assemble(int ndim,
                                        const ConductivityValues & v,
                                        const arma::vec3 & f,
                                        const arma::vec3 & s,
                                        const arma::vec3 & n) const
{
  const double sigma_l = v.sigma_l;
  const double sigma_t = v.sigma_t;
  const double sigma_n = v.sigma_n;
  arma::mat33 I = arma::eye(3, 3);
  arma::mat33 sigma = I;

  if (type_ == ISOTROPIC)
  {
    // diagonal filled with sigma_l (matches original: isotropic uses sigma_l)
    int nr = sigma.n_rows;
    for (int i = 0; i < nr; i++)
      sigma(i, i) = sigma_l;
  }
  else if (type_ == TRANSVERSELY_ISOTROPIC)
  {
    arma::mat tmp = sigma_t * I + (sigma_l - sigma_t) * (f * f.t());
    for (int i = 0; i < ndim; i++)
      for (int j = 0; j < ndim; j++)
        sigma(i, j) = tmp(i, j);
  }
  else if (type_ == ORTHOTROPIC)
  {
    for (int k = 0; k < ndim; k++)
      for (int i = 0; i < ndim; i++)
        sigma(k, i) = sigma_l * f[k] * f[i]
                    + sigma_t * s[k] * s[i]
                    + sigma_n * n[k] * n[i];
  }

  return sigma;
}

PropertyType ConductivityModel::from_string(const std::string & s)
{
  static const std::unordered_map<std::string, PropertyType> lookup = {
    {"isotropic",              ISOTROPIC},
    {"transversely_isotropic", TRANSVERSELY_ISOTROPIC},
    {"orthotropic",            ORTHOTROPIC},
  };

  auto it = lookup.find(s);
  if (it == lookup.end())
    throw std::runtime_error("Unknown conductivity property type: '" + s + "'");
  return it->second;
}

std::string ConductivityModel::to_string(PropertyType t)
{
  switch (t)
  {
    case ISOTROPIC:              return "Isotropic";
    case TRANSVERSELY_ISOTROPIC: return "Transversely_isotropic";
    case ORTHOTROPIC:            return "Orthotropic";
  }
  return "unknown";
}