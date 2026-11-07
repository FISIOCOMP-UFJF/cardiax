#pragma once
/**
 * @file fibrosis.hpp
 * @brief Fibrosis as a change in tissue conductivity, for the monodomain solver.
 *
 * WHAT IT DOES
 * ------------
 * Reads [[fibrosis.region]] blocks from the TOML input and turns them into one
 * conductivity factor per mesh element (a FibrosisField). The solver applies
 * those factors when it builds each element's conductivity tensor.
 *
 * Only the conductivity changes. A node surrounded only by fibrotic elements receives
 * no current from its neighbours, so a ventricular model just stays at rest.
 *
 * Two kinds of remodelling are supported, and a region can use both:
 *
 *   - Diffuse fibrosis: the tissue still conducts, but worse. The fibre,
 *     sheet and normal conductivities are multiplied by `sigma_scale`.
 *
 *   - Patchy fibrosis: a random fraction `density` of the elements becomes
 *     fibrotic. Their conductivity is multiplied by `fibrotic_scale`
 *     (normally 0, so they block conduction).
 *
 * TOML FORMAT
 * -----------
 * Coordinates and lengths are in mesh units (µm in Cardiax).
 *
 *   [fibrosis]
 *   seed = 1508201274              # base seed for regions without their own
 *
 *   [[fibrosis.region]]
 *   name           = "infarct"     # label for the log and error messages
 *   shape          = "box"         # "all" | "sphere" | "box"
 *   min            = [0.0, 1500.0, 7000.0]       # box corners
 *   max            = [3000.0, 5500.0, 13000.0]
 *   # center = [x, y, z]  radius = r              # for a sphere instead
 *   border_zone    = 500.0            # width of the linear fade (default 0)
 *   sigma_scale    = [0.5, 0.3, 0.3]  # fibre, sheet, normal (default 1,1,1)
 *   density        = 0.4              # P(element fibrotic) in the core (default 0)
 *   fibrotic_scale = 0.0              # sigma factor if fibrotic (default 0)
 *   seed           = 42               # optional, per-region pattern
 *
 * Rules:
 *   - Inside the shape the effect is full. Across `border_zone` it fades
 *     linearly to nothing (both the sigma scaling and the fibrosis
 *     probability).
 *   - Regions are applied in the order listed. Overlapping regions multiply
 *     their sigma factors; an element is fibrotic if any region marks it.
 *   - `sigma_scale` acts in the material frame (fibre, sheet, normal), not
 *     along x, y, z. A transversely isotropic conductivity uses only the
 *     first two factors (keep [1] == [2]); an isotropic one only the first.
 *
 * USAGE IN THE SOLVER
 * -------------------
 *   // set_parameters(cfg):
 *   fib_regions = cardiax::read_fibrosis_regions(cfg);
 *
 *   // setup, before the FEM matrices are assembled:
 *   std::vector<arma::vec3> centroids = ...;          // one per element
 *   fib = cardiax::compute_fibrosis_field(fib_regions, centroids);
 *   cardiax::print_fibrosis_summary(fib_regions, fib);
 *
 *   // conductivity tensor of element e:
 *   if (fib.empty()) sigma = cond.at(e, ndim, f, s, n);
 *   else             sigma = cond.at_scaled(e, ndim, f, s, n, fib.sigma_scale[e]);
 *
 * LIMITATIONS
 * -----------
 *   - The random pattern assumes every MPI rank holds the whole mesh
 *     (elements visited in the same global order). With a distributed mesh,
 *     the pattern would depend on the partition.
 *   - Regions are geometric only. Regions tagged in the mesh (e.g. from
 *     imaging) are planned as an extra shape, "marker".
 */

#include <string>
#include <vector>
#include <armadillo>
#include <util/toml.hpp>

namespace cardiax {

/// Geometric shape of a region.
enum class RegionShape
{
  All,     ///< the whole mesh (no geometry keys needed)
  Sphere,  ///< `center` + `radius`
  Box      ///< axis-aligned box from `min` to `max`
};

/**
 * @brief One [[fibrosis.region]] from the TOML file: where it is and what it does.
 *
 * Only the geometry fields of its own shape are used; the others keep their
 * defaults. All lengths are in mesh units.
 */
struct FibrosisRegion
{
  std::string name;   ///< label used in the log and in error messages

  // --- where -------------------------------------------------------------
  RegionShape shape   = RegionShape::All;
  arma::vec3  center  = {0.0, 0.0, 0.0};   ///< sphere center
  double      radius  = 0.0;               ///< sphere radius (> 0)
  arma::vec3  box_min = {0.0, 0.0, 0.0};   ///< box corner with the smallest x, y, z
  arma::vec3  box_max = {0.0, 0.0, 0.0};   ///< box corner with the largest x, y, z
  double      border_zone = 0.0;           ///< width of the linear fade outside the shape

  // --- what --------------------------------------------------------------
  arma::vec3  sigma_scale    = {1.0, 1.0, 1.0};  ///< factors for fibre, sheet, normal conductivity
  double      density        = 0.0;              ///< probability (0..1) that a core element is fibrotic
  double      fibrotic_scale = 0.0;              ///< conductivity factor of a fibrotic element
  unsigned    seed           = 0;                ///< seed of this region's random pattern

  /**
   * @brief How strongly point x is affected by this region.
   * @param x  point in mesh coordinates (normally an element centroid)
   * @return 1 inside the shape, decreasing linearly to 0 across the border
   *         zone, and 0 beyond it. With `border_zone = 0` the result is
   *         exactly 1 or 0.
   */
  double weight(const arma::vec3& x) const;
};

/**
 * @brief Element counts for one region, used to check the result.
 *
 * Compare `fibrotic_in_core / core` with the requested `density`. They differ
 * only by random noise, which is larger when `core` is small.
 */
struct RegionStats
{
  long touched          = 0;   ///< elements with weight > 0 (core + border zone)
  long core             = 0;   ///< elements fully inside the shape (weight = 1)
  long fibrotic_in_core = 0;   ///< core elements this region made fibrotic
};

/**
 * @brief Result of applying all regions to the mesh.
 *
 * Indexed by the element index used when the centroids were passed in, so
 * `sigma_scale[e]` belongs to element e of the mesh.
 */
struct FibrosisField
{
  /// Per element: factors for the fibre, sheet and normal conductivity.
  /// 1 = healthy, 0 = no conduction.
  std::vector<arma::vec3>  sigma_scale;

  /// Per element: true if any region made it fibrotic.
  std::vector<bool>        fibrotic;

  /// Per region, in the same order as the regions vector.
  std::vector<RegionStats> stats;

  /// True when no fibrosis was computed (no [fibrosis] block): the solver
  /// then uses the original conductivity path.
  bool empty() const { return sigma_scale.empty(); }
};

/**
 * @brief Read [fibrosis] and every [[fibrosis.region]] from the TOML input.
 * @param cfg  the parsed TOML file
 * @return the regions in the order they appear; empty if there is no
 *         [fibrosis] block or it has no regions
 * @throws std::runtime_error with the region name and the problem, e.g.
 *         "[fibrosis.region 'x'] 'center' must be a list of 3 numbers"
 */
std::vector<FibrosisRegion> read_fibrosis_regions(const toml::table& cfg);

/**
 * @brief Apply all regions to the mesh.
 * @param regions    as returned by read_fibrosis_regions()
 * @param centroids  centroid of every element, in mesh element order. Every
 *                   MPI rank must pass the same list for the random pattern
 *                   to be identical on all ranks.
 * @return per-element conductivity factors and fibrotic flags, plus
 *         per-region statistics
 */
FibrosisField compute_fibrosis_field(const std::vector<FibrosisRegion>& regions,
                                     const std::vector<arma::vec3>& centroids);

/**
 * @brief Print one line per region (on rank 0 via PetscPrintf).
 *
 * Example:
 *   infarct   1597 elements touched, 697 in core, fibrotic in core 36.2% (requested 40.0%)
 *
 * "0 elements touched" means the region lies outside the mesh: check its
 * coordinates and units against the mesh extent.
 */
void print_fibrosis_summary(const std::vector<FibrosisRegion>& regions,
                            const FibrosisField& field);

} // namespace cardiax
