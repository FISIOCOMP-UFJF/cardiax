/**
 * @file fibrosis.cpp
 * @brief Implementation of the fibrosis module. See fibrosis.hpp for the
 *        TOML format and how the solver uses it.
 *
 * The file has three parts:
 *   1. Geometry   - how far a point is from a region, and its weight.
 *   2. TOML input - parsing and validating [[fibrosis.region]] blocks.
 *   3. The field  - applying the regions to the element centroids.
 */

#include "fibrosis.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <random>
#include <stdexcept>

#include <petsc.h>   // PetscPrintf only (prints once, on rank 0)

namespace cardiax {

// ===========================================================================
//  1. Geometry
// ===========================================================================

/**
 * Distance from point x to the region's shape, in mesh units.
 * Returns 0 for any point inside (or on the surface of) the shape.
 *
 *   All    : every point is inside, so always 0.
 *   Sphere : distance to the center minus the radius, never below 0.
 *   Box    : per axis, how far x lies below `min` or above `max`
 *            (0 if within the box along that axis); the distance is the
 *            length of that "overshoot" vector. This is the exact distance
 *            to the box, including near its edges and corners.
 */
static double distance_to_shape(const FibrosisRegion& region, const arma::vec3& x)
{
  switch (region.shape)
  {
    case RegionShape::All:
      return 0.0;

    case RegionShape::Sphere:
    {
      double dist2 = 0.0;
      for (int k = 0; k < 3; ++k)
        dist2 += (x(k) - region.center(k)) * (x(k) - region.center(k));

      return std::max(0.0, std::sqrt(dist2) - region.radius);
    }

    case RegionShape::Box:
    {
      double dist2 = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        double out = 0.0;   // how far x sticks out of the box along axis k
        if (x(k) < region.box_min(k)) out = region.box_min(k) - x(k);
        if (x(k) > region.box_max(k)) out = x(k) - region.box_max(k);
        dist2 += out * out;
      }
      return std::sqrt(dist2);
    }
  }
  return 0.0;
}

/**
 * Weight w(x) of the region at point x:
 *
 *   w = 1                 inside the shape          (d = 0)
 *   w = 1 - d/border_zone in the border zone         (0 < d < border_zone)
 *   w = 0                 beyond the border zone    (d >= border_zone)
 *
 * With border_zone = 0, the second line never applies and the region has a
 * sharp edge. The weight later scales both the conductivity change and the
 * probability of an element being fibrotic, which gives a smooth transition
 * from damaged to healthy tissue.
 */
double FibrosisRegion::weight(const arma::vec3& x) const
{
  const double d = distance_to_shape(*this, x);

  if (d == 0.0)          return 1.0;   // inside the shape
  if (d >= border_zone)  return 0.0;   // beyond the border zone (or none)
  return 1.0 - d / border_zone;        // fading across the border zone
}

// ===========================================================================
//  2. Reading the TOML file
// ===========================================================================

/**
 * Stop with an error message that names the offending region, e.g.
 *   [fibrosis.region 'infarct'] a sphere needs 'radius' > 0
 * The caller (Monodomain::set_parameters) catches it, prints it and exits.
 */
static void error(const std::string& region, const std::string& message)
{
  throw std::runtime_error("[fibrosis.region '" + region + "'] " + message);
}

/**
 * Read a 3-vector key such as   center = [1.0, 2.0, 3.0]
 * Integers are accepted too (toml++ converts [1, 2, 3] to doubles).
 * Errors if the key is missing, has the wrong length or holds non-numbers.
 */
static arma::vec3 read_vec3(const toml::table& t, const std::string& key,
                            const std::string& region)
{
  const toml::array* a = t[key].as_array();
  if (a == nullptr || a->size() != 3)
    error(region, "'" + key + "' must be a list of 3 numbers");

  arma::vec3 v;
  for (int k = 0; k < 3; ++k)
  {
    std::optional<double> value = (*a)[k].value<double>();
    if (!value)
      error(region, "'" + key + "' must contain only numbers");
    v(k) = *value;
  }
  return v;
}

/**
 * Read the "where" part of a region: `shape`, the geometry keys that shape
 * needs, and `border_zone`. Geometry keys of other shapes are ignored.
 */
static void read_shape(const toml::table& t, FibrosisRegion& r)
{
  const std::string shape = t["shape"].value_or<std::string>("all");

  if (shape == "all")
  {
    r.shape = RegionShape::All;
  }
  else if (shape == "sphere")
  {
    r.shape  = RegionShape::Sphere;
    r.center = read_vec3(t, "center", r.name);
    r.radius = t["radius"].value_or(0.0);
    if (r.radius <= 0.0)
      error(r.name, "a sphere needs 'radius' > 0");
  }
  else if (shape == "box")
  {
    r.shape   = RegionShape::Box;
    r.box_min = read_vec3(t, "min", r.name);
    r.box_max = read_vec3(t, "max", r.name);
    for (int k = 0; k < 3; ++k)
      if (r.box_min(k) > r.box_max(k))
        error(r.name, "box 'min' is larger than 'max'");
  }
  else
  {
    error(r.name, "unknown shape '" + shape + "' (use all, sphere or box)");
  }

  r.border_zone = t["border_zone"].value_or(0.0);
  if (r.border_zone < 0.0)
    error(r.name, "'border_zone' must be >= 0");
}

/**
 * Read the "what" part of a region: `sigma_scale`, `density` and
 * `fibrotic_scale`. Missing keys keep their defaults (no change, no
 * fibrosis, fibrotic elements fully blocked).
 */
static void read_effect(const toml::table& t, FibrosisRegion& r)
{
  if (t.contains("sigma_scale"))
    r.sigma_scale = read_vec3(t, "sigma_scale", r.name);
  for (int k = 0; k < 3; ++k)
    if (r.sigma_scale(k) < 0.0)
      error(r.name, "'sigma_scale' values must be >= 0");

  r.density = t["density"].value_or(0.0);
  if (r.density < 0.0 || r.density > 1.0)
    error(r.name, "'density' must be between 0 and 1");

  r.fibrotic_scale = t["fibrotic_scale"].value_or(0.0);
  if (r.fibrotic_scale < 0.0)
    error(r.name, "'fibrotic_scale' must be >= 0");
}

/**
 * Read every [[fibrosis.region]] in order.
 *
 * Seeds: a region's own `seed` wins. Otherwise region i (counting from 0)
 * gets  [fibrosis].seed + i + 1,  so two regions never share a random pattern
 * even when neither sets a seed. Changing the base seed changes every
 * region's pattern at once.
 */
std::vector<FibrosisRegion> read_fibrosis_regions(const toml::table& cfg)
{
  std::vector<FibrosisRegion> regions;

  const toml::table* fibrosis = cfg["fibrosis"].as_table();
  if (fibrosis == nullptr)
    return regions;                                    // no fibrosis

  const unsigned global_seed = (*fibrosis)["seed"].value_or(0u);

  const toml::array* list = (*fibrosis)["region"].as_array();
  if (list == nullptr)
    return regions;                                    // [fibrosis] without regions

  for (std::size_t i = 0; i < list->size(); ++i)
  {
    const toml::table* t = (*list)[i].as_table();
    if (t == nullptr)
      throw std::runtime_error("each [[fibrosis.region]] must be a table");

    FibrosisRegion r;
    r.name = (*t)["name"].value_or("region_" + std::to_string(i + 1));

    read_shape(*t, r);
    read_effect(*t, r);

    r.seed = (*t)["seed"].value_or(global_seed + static_cast<unsigned>(i) + 1);

    regions.push_back(r);
  }
  return regions;
}

// ===========================================================================
//  3. Applying the regions to the mesh
// ===========================================================================

/**
 * Uniform random number in [0, 1).
 *
 * std::mt19937 produces the same sequence on every compiler (the C++
 * standard fixes it), but std::uniform_real_distribution does not, so the
 * conversion to [0, 1) is written out here. This keeps a given seed giving
 * the same fibrosis pattern with GCC, Clang or any other compiler.
 */
static double random_uniform(std::mt19937& rng)
{
  return rng() / (static_cast<double>(std::mt19937::max()) + 1.0);
}

/**
 * For each region, and for each element e with weight w > 0 at its centroid:
 *
 *   1. Diffuse fibrosis. Each factor k (fibre, sheet, normal) is multiplied by
 *          1 + w * (sigma_scale[k] - 1)
 *      i.e. sigma_scale[k] in the core, 1 (no change) at the outer edge of
 *      the border zone, and in between across it.
 *
 *   2. Patchy fibrosis. The element becomes fibrotic with probability
 *          w * density
 *      and its factors are then multiplied by fibrotic_scale.
 *
 * Factors multiply across regions, so overlapping regions compound
 * (e.g. 0.5 from one region and 0.8 from another gives 0.4). Once a region
 * makes an element fibrotic, it stays fibrotic.
 *
 * Reproducibility: each region has its own generator, seeded with its seed,
 * and draws exactly one number per element, in element order, whether or
 * not the element lies in the region. Element e therefore always receives
 * the same random number for a given seed. Consequences:
 *   - same seed and mesh -> same pattern on every run and every MPI rank
 *     (all ranks hold the whole mesh and visit elements in the same order);
 *   - enlarging a region or raising its density only ADDS fibrotic elements;
 *     the ones already fibrotic stay fibrotic. Useful for parameter sweeps.
 */
FibrosisField compute_fibrosis_field(const std::vector<FibrosisRegion>& regions,
                                     const std::vector<arma::vec3>& centroids)
{
  const std::size_t n_elements = centroids.size();

  // Start from healthy tissue: all factors 1, nothing fibrotic.
  FibrosisField field;
  field.sigma_scale.assign(n_elements, arma::vec3{1.0, 1.0, 1.0});
  field.fibrotic.assign(n_elements, false);
  field.stats.assign(regions.size(), RegionStats{});

  for (std::size_t r = 0; r < regions.size(); ++r)
  {
    const FibrosisRegion& region = regions[r];
    RegionStats&          stats  = field.stats[r];

    std::mt19937 rng(region.seed);

    for (std::size_t e = 0; e < n_elements; ++e)
    {
      // Draw before the weight check: one number per element, always
      // (see "Reproducibility" above).
      const double random = random_uniform(rng);

      const double w = region.weight(centroids[e]);
      if (w == 0.0)
        continue;                                      // element not in this region

      stats.touched++;
      if (w == 1.0)
        stats.core++;

      // 1. Diffuse fibrosis: scale the conductivities, faded by w.
      for (int k = 0; k < 3; ++k)
        field.sigma_scale[e](k) *= 1.0 + w * (region.sigma_scale(k) - 1.0);

      // 2. Patchy fibrosis: fibrotic with probability w * density.
      if (random < w * region.density)
      {
        field.fibrotic[e] = true;
        for (int k = 0; k < 3; ++k)
          field.sigma_scale[e](k) *= region.fibrotic_scale;

        if (w == 1.0)
          stats.fibrotic_in_core++;
      }
    }
  }
  return field;
}

/**
 * One line per region. The fibrotic percentage is shown only for regions
 * with density > 0 and at least one core element. It counts only the core
 * (weight = 1), where the probability is exactly `density`; in the border
 * zone the probability is lower by design.
 */
void print_fibrosis_summary(const std::vector<FibrosisRegion>& regions,
                            const FibrosisField& field)
{
  PetscPrintf(PETSC_COMM_WORLD, " Fibrosis:\n");

  for (std::size_t r = 0; r < regions.size(); ++r)
  {
    const FibrosisRegion& region = regions[r];
    const RegionStats&    stats  = field.stats[r];

    PetscPrintf(PETSC_COMM_WORLD, "   %-20s %ld elements touched, %ld in core",
                region.name.c_str(), stats.touched, stats.core);

    if (region.density > 0.0 && stats.core > 0)
    {
      const double achieved = double(stats.fibrotic_in_core) / stats.core;
      PetscPrintf(PETSC_COMM_WORLD, ", fibrotic in core %.1f%% (requested %.1f%%)",
                  100.0 * achieved, 100.0 * region.density);
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
  }
}

} // namespace cardiax
