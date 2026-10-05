#ifndef KDTREE_UVC_HPP
#define KDTREE_UVC_HPP

#include <string>
#include <vector>

#include "reader_hdf5.hpp"

// =============================================================================
//  Field transfer between cardiac meshes using UVC
//  (Universal Ventricular Coordinates)
// -----------------------------------------------------------------------------
//  The transfer is split in two functions:
//
//    kdtree_uvc_build()     EXPENSIVE, called ONCE. Builds the cylindrical
//                           embedding E = [ab, tm, ab*cos(rt), ab*sin(rt)],
//                           one KD-tree per ventricle and, for every TARGET
//                           node, the list of SOURCE neighbours with their
//                           normalized weights.
//
//    kdtree_uvc_transfer()  CHEAP, called in a LOOP. Only walks the stored
//                           lists: convex mean (continuous) or weighted vote
//                           (categorical). No search, no tree allocation.
//
//  The split works because neighbours and weights depend ONLY on the
//  geometry and the UVC -- never on the field values. Transferring M fields
//  (or M time steps) costs a single k-NN search, done in the build.
//
//  Typical use (what cardiax would do inside the time loop):
//
//      UVCDataTransfer src, tgt;
//      ReaderHDF5      reader_src, reader_tgt;
//      load_uvc_mesh(reader_src, "Patient_1", "source", src);
//      load_uvc_mesh(reader_tgt, "Patient_2", "target", tgt);
//
//      UVCParameters p;             // k=12, gaussian weights, ab clamp on
//      KdtreeUVC kd;
//      kdtree_uvc_build(src, tgt, p, kd);                      // once
//
//      std::vector<double> src_field, tgt_field;
//      for (int step = 0; step < n_steps; step++)
//      {
//        reader_src.read_field_step("vm", step, src_field);
//        kdtree_uvc_transfer(kd, src_field, false, tgt_field);  // in loop
//      }
//
//  Dependencies: C++ standard library, ReaderHDF5 (src/mesh, HighFive) and
//  kdtree-cpp (cdalitz, third party, unmodified).
// =============================================================================

// =============================================================================
//  Mesh + UVC
// =============================================================================

//! Geometry and UVC of one mesh, one value per node. rt ALWAYS in radians.
struct UVCDataTransfer
{
  std::string label;               //!< "source" / "target" (for messages)

  int n_points;
  int n_elements;
  int nen;                         //!< nodes per element

  std::vector<double> xyz;         //!< 3*n_points
  std::vector<int>    tets;        //!< nen*n_elements (connectivity)
  std::vector<double> ab, tm, rt, tv;

  UVCDataTransfer() : label(), n_points(0), n_elements(0), nen(0),
                      xyz(), tets(), ab(), tm(), rt(), tv() {}
};

// =============================================================================
//  Transfer parameters
// =============================================================================

struct UVCParameters
{
  int         k;                   //!< neighbours per target node (12)
  std::string weight;              //!< "gauss" or "idw"
  double      idw_power;           //!< IDW exponent (2.0)
  double      w_ab, w_tm, w_rt;    //!< embedding axis weights (1.0)

  bool        has_tv_split;        //!< false: midpoint of the tv range
  double      tv_split;
  bool        no_clamp_ab;         //!< true: do not clamp target ab

  bool        quiet;               //!< true: build prints nothing

  UVCParameters()
    : k(12), weight("gauss"), idw_power(2.0), w_ab(1.0), w_tm(1.0), w_rt(1.0),
      has_tv_split(false), tv_split(0.0), no_clamp_ab(false),
      quiet(false) {}
};

// =============================================================================
//  Structure built once and queried in a loop
// =============================================================================

//! Everything kdtree_uvc_transfer() needs: the neighbours of each target node
//! and their weights. No KD-tree survives the build -- that is why the
//! struct is copyable and cheap to keep around.
//!
//! The lists are stored in fixed-width rows of size k_max:
//!   neighbors[k_max*i + j] -> SOURCE node (valid only for j < n_neighbors[i])
//!   weights  [k_max*i + j] -> normalized weight (each row sums to 1)
struct KdtreeUVC
{
  int n_source;                    //!< nodes of the source mesh
  int n_target;                    //!< nodes of the target mesh
  int k_max;                       //!< row width of neighbors/weights

  std::vector<int>    neighbors;   //!< k_max*n_target
  std::vector<double> weights;     //!< k_max*n_target
  std::vector<int>    n_neighbors; //!< n_target

  //! Target nodes that got no neighbour at all (non-finite UVC, or a
  //! ventricle with no source nodes). They copy the value of another target
  //! node: target[copy_to[t]] = target[copy_from[t]], where copy_from is the
  //! closest node that did get neighbours.
  std::vector<int> copy_to;
  std::vector<int> copy_from;

  double tv_split;                 //!< LV/RV threshold actually used
  int    n_unfilled;               //!< target nodes left as NaN anyway

  KdtreeUVC()
    : n_source(0), n_target(0), k_max(0), neighbors(), weights(),
      n_neighbors(), copy_to(), copy_from(), tv_split(0.5), n_unfilled(0) {}
};

// =============================================================================
//  Reading
// =============================================================================

//! Reads geometry + UVC from an HDF5 file and leaves rt in radians.
//! The reader is left OPEN on purpose: the field to transfer usually comes
//! from the same file and is read later, step by step.
//! Returns false (message in cerr) if any UVC field is missing.
bool load_uvc_mesh(ReaderHDF5 & reader, const std::string & filename,
                   const std::string & label, UVCDataTransfer & m,
                   const std::string & n_ab = "ab",
                   const std::string & n_tm = "tm",
                   const std::string & n_rt = "rt",
                   const std::string & n_tv = "tv");

//! Reads a scalar field as stored in the file (per node or per element).
//! name: just the name ("tecido") or the path ("vertex_field/tecido").
//! cell_data and path are filled on return; path is the dataset PATH in the
//! file ("/vertex_field/tecido"), unique even if the name exists in both
//! groups.
bool read_field_uvc(ReaderHDF5 & reader, const std::string & name,
                    std::vector<double> & values, bool & cell_data,
                    std::string & path);

// =============================================================================
//  Main interface
// =============================================================================

//! CALLED ONCE. Builds the embedding, the trees and the neighbour lists.
//! Looks at no field: only xyz/ab/tm/rt/tv of the two meshes.
//! Returns false if the meshes are inconsistent or empty.
bool kdtree_uvc_build(const UVCDataTransfer & source,
                      const UVCDataTransfer & target,
                      const UVCParameters & p, KdtreeUVC & kd);

//! CALLED IN A LOOP. source_field has one value per SOURCE node;
//! target_field gets one value per TARGET node (NaN where undecidable).
//!
//!   categorical = false -> convex mean of the neighbours (no overshoot)
//!   categorical = true  -> weighted vote, preserves 0/1 and labels
//!
//! Neighbours with non-finite values are skipped and the remaining weights
//! renormalized -- which is why the build does not need to know the field.
bool kdtree_uvc_transfer(const KdtreeUVC & kd,
                         const std::vector<double> & source_field,
                         bool categorical,
                         std::vector<double> & target_field);

// =============================================================================
//  PointData <-> CellData conversion (keeps the field type)
// =============================================================================

//! Distinct finite values, in increasing order
std::vector<double> distinct_values_uvc(const std::vector<double> & v);

//! Field per CELL -> per NODE.
//! continuous  -> mean of the incident cells
//! categorical -> majority vote (binary: fraction > 0.5; multiclass: mode)
std::vector<double> cell_to_node_uvc(const UVCDataTransfer & m,
                                     const std::vector<double> & cell_values,
                                     bool categorical);

//! Field per NODE -> per CELL.
//! continuous  -> mean of the cell's nodes
//! categorical -> binary: marks the cell when the fraction of nodes in the
//!                positive class reaches frac; multiclass: mode of the nodes
std::vector<double> node_to_cell_uvc(const UVCDataTransfer & m,
                                     const std::vector<double> & node_values,
                                     bool categorical, double frac = 0.5);

// =============================================================================
//  Output
// =============================================================================

struct OutputField
{
  std::string name;
  std::vector<double> values;
  bool cell_data;                  //!< false = PointData, true = CellData

  OutputField() : name(), values(), cell_data(false) {}
};

//! Writes the mesh and the fields to an ASCII .vtu, hand-written (no VTK).
//! Each field goes to PointData or CellData according to cell_data.
bool save_vtu_uvc(const std::string & filename, const UVCDataTransfer & m,
                  const std::vector<OutputField> & fields);

#endif /* KDTREE_UVC_HPP */
