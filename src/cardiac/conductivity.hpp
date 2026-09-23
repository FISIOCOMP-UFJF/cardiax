#ifndef CONDUCTIVITY_MODEL_HPP
#define CONDUCTIVITY_MODEL_HPP

#include <map>
#include <string>
#include <armadillo>

// Anisotropy type for the conductivity tensor. Owned here rather than in the
// mesh, since the material model (not the geometry) defines it.
enum PropertyType
{
  ISOTROPIC,
  TRANSVERSELY_ISOTROPIC,
  ORTHOTROPIC
};

// Conductivity magnitudes in the three material directions
// (longitudinal/fiber, transverse/sheet, normal).
struct ConductivityValues
{
  double sigma_l = 0.0;
  double sigma_t = 0.0;
  double sigma_n = 0.0;
};

class ConductivityModel
{
public:
  ConductivityModel() = default;

  void set_type(PropertyType t) { type_ = t; }
  PropertyType type() const { return type_; }
  void set_default(const ConductivityValues & v) { default_ = v; }
  void set_region(int region_id, const ConductivityValues & v) { table_[region_id] = v; }
  const ConductivityValues & values(int region_id) const;

  arma::mat33 at(int region_id, int ndim,
                 const arma::vec3 & f,
                 const arma::vec3 & s,
                 const arma::vec3 & n) const;

  arma::mat33 assemble(int ndim,
                       const ConductivityValues & v,
                       const arma::vec3 & f,
                       const arma::vec3 & s,
                       const arma::vec3 & n) const;

  static PropertyType from_string(const std::string & s);
  static std::string  to_string(PropertyType t);

private:
  PropertyType                      type_ = ISOTROPIC;
  ConductivityValues                default_;
  std::map<int, ConductivityValues> table_;
};

#endif // CONDUCTIVITY_MODEL_HPP