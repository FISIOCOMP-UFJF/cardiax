#include "kdtree_uvc.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "kdtree.hpp"      // cdalitz/kdtree-cpp, unmodified

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// =============================================================================
//  Local helpers (anonymous namespace: they do not leak to the rest of the
//  project)
// =============================================================================

namespace
{
  const double PI_     = 3.1415926535897932384626433832795;
  const double TWO_PI  = 2.0 * PI_;
  const double NOT_NUM = 0.0 / 0.0;

  //! true if the value is neither NaN nor infinite
  bool is_finite(double v)
  {
    return (v == v) && (v > -1e300) && (v < 1e300);
  }

  //! Cylindrical embedding of one node: E = [ab*w_ab, tm*w_tm,
  //!                                         ab*cos(rt)*w_rt, ab*sin(rt)*w_rt]
  //! cos/sin handle the +-pi seam without ghost copies, and the radius
  //! ab -> 0 makes rt stop mattering at the apex, where it degenerates.
  void embed(double ab, double tm, double rt, const UVCParameters & p,
             double e[4])
  {
    double cx = ab * std::cos(rt);
    double cy = ab * std::sin(rt);
    if (!is_finite(cx) || !is_finite(cy)) { cx = 0.0; cy = 0.0; }  // rt NaN

    e[0] = ab * p.w_ab;
    e[1] = tm * p.w_tm;
    e[2] = cx * p.w_rt;
    e[3] = cy * p.w_rt;
  }

  //! Nearest neighbour of q in a point cloud of dimension dim.
  //! Returns the position of the point in the list, or -1.
  int nearest(Kdtree::KdTree & tree, const double * q, int dim)
  {
    Kdtree::KdNodeVector nb;
    tree.k_nearest_neighbors(Kdtree::CoordPoint(q, q + dim), 1, &nb);
    if (nb.empty()) return -1;
    return nb[0].index;
  }

}  // anonymous namespace

// =============================================================================
//  Reading
// =============================================================================

bool read_field_uvc(ReaderHDF5 & reader, const string & name,
                    vector<double> & values, bool & cell_data,
                    string & path)
{
  // name may be just the name ("tecido") or the path ("vertex_field/tecido")
  const int idx = reader.find_field(name);
  if (idx < 0) return false;

  const FieldInfo & f = reader.get_field(idx);
  if (f.n_comp != 1)
  {
    cerr << "[ERROR] field '" << f.name << "' has " << f.n_comp
         << " components; expected a scalar." << endl;
    return false;
  }

  if (!reader.read_field_step(f.path, 0, values)) return false;

  cell_data = f.cell_centered;
  path      = f.path;           // path: unique key for later reads
  return true;
}

// -----------------------------------------------------------------------------

bool load_uvc_mesh(ReaderHDF5 & reader, const string & filename,
                   const string & label, UVCDataTransfer & m,
                   const string & n_ab, const string & n_tm,
                   const string & n_rt, const string & n_tv)
{
  if (!reader.open(filename))
  {
    cerr << "[ERROR] cannot open '" << filename << "'" << endl;
    return false;
  }

  m.label      = label;
  m.n_points   = reader.get_n_points();
  m.n_elements = reader.get_n_elements();
  m.nen        = reader.get_nen();
  m.xyz        = reader.get_coordinates();
  m.tets       = reader.get_connectivity();

  const string names[4] = { n_ab, n_tm, n_rt, n_tv };
  vector<double> * dest[4] = { &m.ab, &m.tm, &m.rt, &m.tv };
  string missing;

  for (int c = 0; c < 4; c++)
  {
    bool cell_data = false;
    string path;
    vector<double> v;

    if (!read_field_uvc(reader, names[c], v, cell_data, path))
    {
      missing += (missing.empty() ? "" : ", ") + names[c];
      continue;
    }

    if (cell_data)
    {
      cout << "  [WARNING] UVC '" << names[c] << "' is per ELEMENT; "
           << "converted to nodes by averaging." << endl;
      v = cell_to_node_uvc(m, v, false);
    }
    *dest[c] = v;
  }

  if (!missing.empty())
  {
    cerr << "[ERROR] UVC fields missing in '" << filename << "': "
         << missing << endl;
    cerr << "        available fields:";
    for (int i = 0; i < reader.get_n_fields(); i++)
      cerr << " " << reader.get_field(i).name;
    cerr << endl;
    return false;
  }

  // ------------------------------------------------------------- rt unit
  // The embedding uses cos(rt)/sin(rt): rt must be in RADIANS.
  double lo = 1e300, hi = -1e300;
  for (int i = 0; i < m.n_points; i++)
  {
    if (!is_finite(m.rt[(size_t) i])) continue;
    if (m.rt[(size_t) i] < lo) lo = m.rt[(size_t) i];
    if (m.rt[(size_t) i] > hi) hi = m.rt[(size_t) i];
  }

  if (!(hi > 1.6 || lo < -1.6))
  {
    cout << "  [INFO] " << label << ": rt seems to be in turns [" << lo
         << ", " << hi << "]; converting to radians (rt * 2*pi)." << endl;
    for (int i = 0; i < m.n_points; i++)
      if (is_finite(m.rt[(size_t) i])) m.rt[(size_t) i] *= TWO_PI;
  }

  cout << "  " << label << " : " << m.n_points << " nodes, "
       << m.n_elements << " cells (" << m.nen << " nodes each)" << endl;

  return true;
}

// =============================================================================
//  Build -- CALLED ONCE
// =============================================================================

bool kdtree_uvc_build(const UVCDataTransfer & source,
                      const UVCDataTransfer & target,
                      const UVCParameters & p, KdtreeUVC & kd)
{
  const int ns = source.n_points;
  const int nt = target.n_points;

  if (ns <= 0 || nt <= 0)
  {
    cerr << "[ERROR] kdtree_uvc_build: empty mesh (source " << ns
         << " nodes, target " << nt << " nodes)." << endl;
    return false;
  }

  if ((int) source.ab.size() != ns || (int) source.tm.size() != ns ||
      (int) source.rt.size() != ns || (int) source.tv.size() != ns ||
      (int) target.ab.size() != nt || (int) target.tm.size() != nt ||
      (int) target.rt.size() != nt || (int) target.tv.size() != nt)
  {
    cerr << "[ERROR] kdtree_uvc_build: UVC size differs from the number "
         << "of nodes." << endl;
    return false;
  }

  if (p.k < 1)
  {
    cerr << "[ERROR] kdtree_uvc_build: k must be >= 1." << endl;
    return false;
  }
  if (p.weight != "gauss" && p.weight != "idw")
  {
    cerr << "[ERROR] kdtree_uvc_build: weight must be 'gauss' or 'idw'."
         << endl;
    return false;
  }

  kd.n_source = ns;
  kd.n_target = nt;

  // --------------------------------------------------------------- 1) masks
  // Note that the field does NOT enter here: validity depends only on the
  // UVC. Neighbours with non-finite values are skipped later, in
  // kdtree_uvc_transfer().
  vector<char> ok_s((size_t) ns, 0);
  vector<char> ok_t((size_t) nt, 0);
  int n_ok_s = 0, n_ok_t = 0;

  for (int i = 0; i < ns; i++)
  {
    ok_s[(size_t) i] = (char) (is_finite(source.ab[(size_t) i]) &&
                               is_finite(source.tm[(size_t) i]) &&
                               is_finite(source.tv[(size_t) i]));
    n_ok_s += ok_s[(size_t) i];
  }
  for (int i = 0; i < nt; i++)
  {
    ok_t[(size_t) i] = (char) (is_finite(target.ab[(size_t) i]) &&
                               is_finite(target.tm[(size_t) i]) &&
                               is_finite(target.tv[(size_t) i]));
    n_ok_t += ok_t[(size_t) i];
  }

  if (n_ok_s == 0)
  {
    cerr << "[ERROR] kdtree_uvc_build: no valid node in the source." << endl;
    return false;
  }

  if (!p.quiet)
    cout << "  valid source nodes: " << n_ok_s << "/" << ns
         << " | valid target nodes: " << n_ok_t << "/" << nt << endl;

  // ----------------------------------------------------------- 2) tv_split
  if (p.has_tv_split) kd.tv_split = p.tv_split;
  else
  {
    double lo = 1e300, hi = -1e300;
    for (int i = 0; i < ns; i++)
      if (is_finite(source.tv[(size_t) i]))
      {
        if (source.tv[(size_t) i] < lo) lo = source.tv[(size_t) i];
        if (source.tv[(size_t) i] > hi) hi = source.tv[(size_t) i];
      }
    for (int i = 0; i < nt; i++)
      if (is_finite(target.tv[(size_t) i]))
      {
        if (target.tv[(size_t) i] < lo) lo = target.tv[(size_t) i];
        if (target.tv[(size_t) i] > hi) hi = target.tv[(size_t) i];
      }
    kd.tv_split = (lo > hi) ? 0.5 : 0.5 * (lo + hi);
  }
  const double tv_split = kd.tv_split;

  if (!p.quiet)
    cout << "  TV_SPLIT = " << tv_split << "  (LV: tv<split, RV: tv>=split)"
         << endl;

  // ------------------------------------------------ 3) clamp of target ab
  vector<double> ab_t = target.ab;
  if (!p.no_clamp_ab)
  {
    double lo = 1e300, hi = -1e300;
    for (int i = 0; i < ns; i++)
      if (ok_s[(size_t) i])
      {
        if (source.ab[(size_t) i] < lo) lo = source.ab[(size_t) i];
        if (source.ab[(size_t) i] > hi) hi = source.ab[(size_t) i];
      }

    for (int i = 0; i < nt; i++)
      if (is_finite(ab_t[(size_t) i]))
      {
        if (ab_t[(size_t) i] < lo) ab_t[(size_t) i] = lo;
        if (ab_t[(size_t) i] > hi) ab_t[(size_t) i] = hi;
      }

    if (!p.quiet)
      cout << "  target ab clamped to [" << lo << ", " << hi << "]." << endl;
  }

  // --------------------------------------------- 4) cylindrical embedding
  vector<double> E_s((size_t) 4 * ns, 0.0);
  vector<double> E_t((size_t) 4 * nt, 0.0);

  for (int i = 0; i < ns; i++)
    embed(source.ab[(size_t) i], source.tm[(size_t) i], source.rt[(size_t) i],
          p, &E_s[(size_t) 4 * i]);

  for (int i = 0; i < nt; i++)
    embed(ab_t[(size_t) i], target.tm[(size_t) i], target.rt[(size_t) i],
          p, &E_t[(size_t) 4 * i]);

  // ------------------------------- 5) number of source nodes per ventricle
  // Must come before allocating the rows: k_max is the largest effective k
  // of the two ventricles.
  int n_src[2] = { 0, 0 };
  for (int i = 0; i < ns; i++)
  {
    if (!ok_s[(size_t) i]) continue;
    n_src[(source.tv[(size_t) i] < tv_split) ? 0 : 1]++;
  }

  int kk[2];
  for (int v = 0; v < 2; v++)
    kk[v] = (p.k < n_src[v]) ? p.k : n_src[v];

  kd.k_max = (kk[0] > kk[1]) ? kk[0] : kk[1];
  if (kd.k_max < 1)
  {
    cerr << "[ERROR] kdtree_uvc_build: no ventricle with source nodes."
         << endl;
    return false;
  }

  kd.neighbors.assign((size_t) kd.k_max * nt, -1);
  kd.weights.assign((size_t) kd.k_max * nt, 0.0);
  kd.n_neighbors.assign((size_t) nt, 0);

  const bool gauss = (p.weight == "gauss");

  // ----------------------------------------------- 6) k-NN per ventricle
  for (int vent = 0; vent < 2; vent++)
  {
    const string vname = vent ? "RV" : "LV";

    // --- source nodes of this ventricle (KdNode::index = position in id_s)
    Kdtree::KdNodeVector pts;
    vector<int> id_s;
    for (int i = 0; i < ns; i++)
    {
      if (!ok_s[(size_t) i]) continue;
      const bool lv = (source.tv[(size_t) i] < tv_split);
      if (lv == (vent == 1)) continue;
      const double * e = &E_s[(size_t) 4 * i];
      pts.push_back(Kdtree::KdNode(Kdtree::CoordPoint(e, e + 4), NULL,
                                   (int) id_s.size()));
      id_s.push_back(i);
    }

    // --- target nodes of this ventricle ---
    vector<int> id_t;
    for (int i = 0; i < nt; i++)
    {
      if (!ok_t[(size_t) i]) continue;
      const bool lv = (target.tv[(size_t) i] < tv_split);
      if (lv == (vent == 1)) continue;
      id_t.push_back(i);
    }

    if (id_s.empty() || id_t.empty())
    {
      if (!p.quiet)
        cout << "  " << vname << ": no valid nodes, skipping." << endl;
      continue;
    }

    // id_s is not empty: Kdtree::KdTree accepts the list
    Kdtree::KdTree tree(&pts);
    Kdtree::KdNodeVector().swap(pts);   // the tree keeps its own copy

    const int k = kk[vent];
    Kdtree::CoordPoint   q(4);
    Kdtree::KdNodeVector nb;            // in increasing distance order
    vector<double> d((size_t) k, 0.0);
    vector<double> w((size_t) k, 0.0);

    for (size_t t = 0; t < id_t.size(); t++)
    {
      const int i = id_t[t];
      for (int c = 0; c < 4; c++) q[(size_t) c] = E_t[(size_t) 4 * i + c];

      tree.k_nearest_neighbors(q, (size_t) k, &nb);
      const int nv = (int) nb.size();
      if (nv == 0) continue;

      // kdtree-cpp does not return distances: recompute them from the point
      for (int j = 0; j < nv; j++)
      {
        const Kdtree::CoordPoint & pt = nb[(size_t) j].point;
        double s2 = 0.0;
        for (int c = 0; c < 4; c++)
        {
          const double dc = q[(size_t) c] - pt[(size_t) c];
          s2 += dc * dc;
        }
        d[(size_t) j] = std::sqrt(s2);
      }

      // --- weights ---
      if (gauss)
      {
        double mean = 0.0;
        for (int j = 0; j < nv; j++) mean += d[(size_t) j];
        mean /= (double) nv;
        const double h = (mean > 1e-9) ? mean : 1e-9;
        for (int j = 0; j < nv; j++)
        {
          const double z = d[(size_t) j] / h;
          w[(size_t) j] = std::exp(-z * z);
        }
      }
      else
      {
        for (int j = 0; j < nv; j++)
        {
          const double dd = (d[(size_t) j] > 1e-12) ? d[(size_t) j] : 1e-12;
          w[(size_t) j] = 1.0 / std::pow(dd, p.idw_power);
        }
      }

      // exact coincidence: only the first neighbour counts
      if (d[0] < 1e-12)
      {
        for (int j = 0; j < nv; j++) w[(size_t) j] = 0.0;
        w[0] = 1.0;
      }

      double sum = 0.0;
      for (int j = 0; j < nv; j++) sum += w[(size_t) j];
      if (sum <= 0.0) { w[0] = 1.0; sum = 1.0; }

      // --- store the row already normalized ---
      const size_t base = (size_t) kd.k_max * i;
      for (int j = 0; j < nv; j++)
      {
        kd.neighbors[base + j] = id_s[(size_t) nb[(size_t) j].index];
        kd.weights[base + j]   = w[(size_t) j] / sum;
      }
      kd.n_neighbors[(size_t) i] = nv;
    }

    if (!p.quiet)
      cout << "  " << vname << ": " << id_t.size() << " target nodes (k="
           << k << ")." << endl;
  }

  // ------------------- 7) target nodes without neighbours: copy from a peer
  // Who copies from whom depends only on the geometry, so it is decided here
  // and applied in every transfer.
  kd.copy_to.clear();
  kd.copy_from.clear();
  kd.n_unfilled = 0;

  {
    Kdtree::KdNodeVector pts_e;       // filled nodes, in embedding space
    Kdtree::KdNodeVector pts_x;       // the same nodes, in physical space
    vector<int> filled;
    vector<int> unfilled;

    for (int i = 0; i < nt; i++)
    {
      if (kd.n_neighbors[(size_t) i] > 0)
      {
        const double * e = &E_t[(size_t) 4 * i];
        const double * x = &target.xyz[(size_t) 3 * i];
        pts_e.push_back(Kdtree::KdNode(Kdtree::CoordPoint(e, e + 4), NULL,
                                       (int) filled.size()));
        pts_x.push_back(Kdtree::KdNode(Kdtree::CoordPoint(x, x + 3), NULL,
                                       (int) filled.size()));
        filled.push_back(i);
      }
      else unfilled.push_back(i);
    }

    if (!unfilled.empty() && !filled.empty())
    {
      Kdtree::KdTree tree_e(&pts_e);

      // The physical tree is built only if some orphan node has a
      // non-finite embedding (ab or tm NaN): querying the embedding tree
      // with NaN would return garbage.
      Kdtree::KdTree * tree_x = 0;

      for (size_t t = 0; t < unfilled.size(); t++)
      {
        const int i = unfilled[t];
        const double * e = &E_t[(size_t) 4 * i];

        bool e_ok = true;
        for (int c = 0; c < 4; c++) if (!is_finite(e[c])) e_ok = false;

        int pos = -1;
        if (e_ok)
        {
          pos = nearest(tree_e, e, 4);
        }
        else
        {
          if (tree_x == 0) tree_x = new Kdtree::KdTree(&pts_x);
          pos = nearest(*tree_x, &target.xyz[(size_t) 3 * i], 3);
        }

        if (pos < 0) { kd.n_unfilled++; continue; }

        kd.copy_to.push_back(i);
        kd.copy_from.push_back(filled[(size_t) pos]);
      }

      delete tree_x;

      if (!p.quiet)
        cout << "  " << kd.copy_to.size()
             << " target nodes filled from a UVC neighbour." << endl;
    }
    else kd.n_unfilled = (int) unfilled.size();
  }

  return true;
}

// =============================================================================
//  Transfer -- CALLED IN A LOOP
// =============================================================================

bool kdtree_uvc_transfer(const KdtreeUVC & kd,
                         const vector<double> & source_field,
                         bool categorical,
                         vector<double> & target_field)
{
  if ((int) source_field.size() != kd.n_source)
  {
    cerr << "[ERROR] kdtree_uvc_transfer: field has " << source_field.size()
         << " values; the source has " << kd.n_source << " nodes." << endl;
    return false;
  }
  if (kd.k_max < 1 || (int) kd.n_neighbors.size() != kd.n_target)
  {
    cerr << "[ERROR] kdtree_uvc_transfer: structure not built "
         << "(call kdtree_uvc_build first)." << endl;
    return false;
  }

  target_field.assign((size_t) kd.n_target, NOT_NUM);

  const int kmax = kd.k_max;

  for (int i = 0; i < kd.n_target; i++)
  {
    const int nv = kd.n_neighbors[(size_t) i];
    if (nv <= 0) continue;

    const size_t base = (size_t) kmax * i;

    // ---------------------------------------------------------- continuous
    if (!categorical)
    {
      double v = 0.0, sum = 0.0;
      for (int j = 0; j < nv; j++)
      {
        const double val = source_field[(size_t) kd.neighbors[base + j]];
        if (!is_finite(val)) continue;          // neighbour without value
        v   += kd.weights[base + j] * val;
        sum += kd.weights[base + j];
      }
      if (sum > 0.0) target_field[(size_t) i] = v / sum;   // renormalize
      continue;
    }

    // --------------------------------------------------------- categorical
    // Weighted vote over the values present among the neighbours. Since nv
    // is small (k ~ 12), scanning for repeats is cheaper than keeping a
    // global class list -- and avoids passing it to this function.
    // Ties go to the SMALLER label.
    double best_w = -1.0;
    double best_v = NOT_NUM;

    for (int j = 0; j < nv; j++)
    {
      const double val = source_field[(size_t) kd.neighbors[base + j]];
      if (!is_finite(val)) continue;

      bool repeated = false;
      for (int j2 = 0; j2 < j && !repeated; j2++)
        if (source_field[(size_t) kd.neighbors[base + j2]] == val)
          repeated = true;
      if (repeated) continue;

      double sw = 0.0;
      for (int j2 = 0; j2 < nv; j2++)
        if (source_field[(size_t) kd.neighbors[base + j2]] == val)
          sw += kd.weights[base + j2];

      if (sw > best_w || (sw == best_w && val < best_v))
      {
        best_w = sw;
        best_v = val;
      }
    }

    if (best_w >= 0.0) target_field[(size_t) i] = best_v;
  }

  // ------------------------------- copies (target nodes without neighbours)
  for (size_t t = 0; t < kd.copy_to.size(); t++)
    target_field[(size_t) kd.copy_to[t]] =
        target_field[(size_t) kd.copy_from[t]];

  return true;
}

// =============================================================================
//  PointData <-> CellData conversion
// =============================================================================

vector<double> distinct_values_uvc(const vector<double> & v)
{
  vector<double> c;
  for (size_t i = 0; i < v.size(); i++)
    if (is_finite(v[i])) c.push_back(v[i]);

  std::sort(c.begin(), c.end());
  c.erase(std::unique(c.begin(), c.end()), c.end());
  return c;
}

// -----------------------------------------------------------------------------

vector<double> cell_to_node_uvc(const UVCDataTransfer & m,
                                const vector<double> & cell_values,
                                bool categorical)
{
  const int np  = m.n_points;
  const int ne  = m.n_elements;
  const int nen = m.nen;

  const vector<double> classes = distinct_values_uvc(cell_values);

  // ---------------------------------------------------------- continuous
  if (!categorical)
  {
    vector<double> sum((size_t) np, 0.0);
    vector<double> cnt((size_t) np, 0.0);

    for (int e = 0; e < ne; e++)
      for (int j = 0; j < nen; j++)
      {
        const int node = m.tets[(size_t) nen * e + j];
        sum[(size_t) node] += cell_values[(size_t) e];
        cnt[(size_t) node] += 1.0;
      }

    for (int i = 0; i < np; i++)
      if (cnt[(size_t) i] > 0.0) sum[(size_t) i] /= cnt[(size_t) i];
    return sum;
  }

  // -------------------------------------------------- categorical binary
  if (classes.size() <= 2)
  {
    const double low  = classes.empty() ? 0.0 : classes.front();
    const double high = (classes.size() == 2) ? classes.back() : low;

    vector<double> sum((size_t) np, 0.0);
    vector<double> cnt((size_t) np, 0.0);

    for (int e = 0; e < ne; e++)
    {
      const double b = (classes.size() == 2 && cell_values[(size_t) e] > low)
                       ? 1.0 : 0.0;
      for (int j = 0; j < nen; j++)
      {
        const int node = m.tets[(size_t) nen * e + j];
        sum[(size_t) node] += b;
        cnt[(size_t) node] += 1.0;
      }
    }

    vector<double> out((size_t) np, low);
    for (int i = 0; i < np; i++)
    {
      const double f = cnt[(size_t) i] ? sum[(size_t) i] / cnt[(size_t) i]
                                       : 0.0;
      // Exact tie (f = 0.5) goes to the SMALLER class, matching numpy's
      // round-half-to-even used by the original Python script.
      out[(size_t) i] = (f > 0.5) ? high : low;
    }
    return out;
  }

  // ---------------------------------------------- categorical multiclass
  const size_t nc = classes.size();
  vector<double> count((size_t) np * nc, 0.0);

  for (int e = 0; e < ne; e++)
  {
    const size_t c = (size_t) (std::lower_bound(classes.begin(), classes.end(),
                                                cell_values[(size_t) e])
                               - classes.begin());
    if (c >= nc) continue;
    for (int j = 0; j < nen; j++)
      count[(size_t) m.tets[(size_t) nen * e + j] * nc + c] += 1.0;
  }

  vector<double> out((size_t) np, classes[0]);
  for (int i = 0; i < np; i++)
  {
    size_t best = 0;
    for (size_t c = 1; c < nc; c++)
      if (count[(size_t) i * nc + c] > count[(size_t) i * nc + best])
        best = c;
    out[(size_t) i] = classes[best];
  }
  return out;
}

// -----------------------------------------------------------------------------

vector<double> node_to_cell_uvc(const UVCDataTransfer & m,
                                const vector<double> & node_values,
                                bool categorical, double frac)
{
  const int ne  = m.n_elements;
  const int nen = m.nen;

  vector<double> out((size_t) ne, 0.0);

  // ---------------------------------------------------------- continuous
  if (!categorical)
  {
    for (int e = 0; e < ne; e++)
    {
      double sum = 0.0;
      int    n   = 0;
      for (int j = 0; j < nen; j++)
      {
        const double v = node_values[(size_t) m.tets[(size_t) nen * e + j]];
        if (!is_finite(v)) continue;
        sum += v;
        n++;
      }
      out[(size_t) e] = n ? sum / (double) n : NOT_NUM;
    }
    return out;
  }

  const vector<double> classes = distinct_values_uvc(node_values);

  // ---- binary with zero as the negative class: fraction of positive nodes
  const bool binary = (classes.size() <= 2 && !classes.empty() &&
                       classes.front() == 0.0);

  if (binary)
  {
    const double positive = (classes.size() == 2) ? classes.back() : 1.0;

    for (int e = 0; e < ne; e++)
    {
      int pos = 0;
      for (int j = 0; j < nen; j++)
        if (node_values[(size_t) m.tets[(size_t) nen * e + j]] == positive)
          pos++;
      const double f = (double) pos / (double) nen;
      out[(size_t) e] = (f >= frac) ? positive : 0.0;
    }
    return out;
  }

  // ------------------------------------------- multiclass: mode of the nodes
  for (int e = 0; e < ne; e++)
  {
    double best   = node_values[(size_t) m.tets[(size_t) nen * e]];
    int    best_n = 0;

    for (int j = 0; j < nen; j++)
    {
      const double v = node_values[(size_t) m.tets[(size_t) nen * e + j]];
      int n = 0;
      for (int j2 = 0; j2 < nen; j2++)
        if (node_values[(size_t) m.tets[(size_t) nen * e + j2]] == v) n++;
      if (n > best_n) { best_n = n; best = v; }
    }
    out[(size_t) e] = best;
  }
  return out;
}

// =============================================================================
//  Writing the target mesh to .vtu (ASCII, hand-written -- no VTK)
// =============================================================================

bool save_vtu_uvc(const string & filename, const UVCDataTransfer & m,
                  const vector<OutputField> & fields)
{
  std::ofstream f(filename.c_str());
  if (!f)
  {
    cerr << "[ERROR] cannot write '" << filename << "'" << endl;
    return false;
  }

  // VTK cell type from the number of nodes per element
  int vtk_type = 10;                       // VTK_TETRA
  if      (m.nen == 3) vtk_type = 5;       // VTK_TRIANGLE
  else if (m.nen == 8) vtk_type = 12;      // VTK_HEXAHEDRON
  else if (m.nen == 2) vtk_type = 3;       // VTK_LINE

  f << std::setprecision(10);
  f << "<?xml version=\"1.0\"?>\n";
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    << "byte_order=\"LittleEndian\">\n";
  f << "  <UnstructuredGrid>\n";
  f << "    <Piece NumberOfPoints=\"" << m.n_points
    << "\" NumberOfCells=\"" << m.n_elements << "\">\n";

  // ------------------------------------------------------------- points
  f << "      <Points>\n";
  f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
    << "format=\"ascii\">\n";
  for (int i = 0; i < m.n_points; i++)
    f << "          " << m.xyz[(size_t) 3 * i + 0] << " "
      << m.xyz[(size_t) 3 * i + 1] << " " << m.xyz[(size_t) 3 * i + 2] << "\n";
  f << "        </DataArray>\n      </Points>\n";

  // -------------------------------------------------------------- cells
  f << "      <Cells>\n";
  f << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
    << "format=\"ascii\">\n";
  for (int e = 0; e < m.n_elements; e++)
  {
    f << "          ";
    for (int j = 0; j < m.nen; j++) f << m.tets[(size_t) m.nen * e + j] << " ";
    f << "\n";
  }
  f << "        </DataArray>\n";

  f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
    << "          ";
  for (int e = 0; e < m.n_elements; e++) f << (m.nen * (e + 1)) << " ";
  f << "\n        </DataArray>\n";

  f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
    << "          ";
  for (int e = 0; e < m.n_elements; e++) f << vtk_type << " ";
  f << "\n        </DataArray>\n      </Cells>\n";

  // ------------------------------------------------------------- fields
  f << "      <PointData>\n";
  for (size_t c = 0; c < fields.size(); c++)
  {
    if (fields[c].cell_data) continue;
    f << "        <DataArray type=\"Float64\" Name=\"" << fields[c].name
      << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < fields[c].values.size(); i++)
      f << fields[c].values[i] << " ";
    f << "\n        </DataArray>\n";
  }
  f << "      </PointData>\n";

  f << "      <CellData>\n";
  for (size_t c = 0; c < fields.size(); c++)
  {
    if (!fields[c].cell_data) continue;
    f << "        <DataArray type=\"Float64\" Name=\"" << fields[c].name
      << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < fields[c].values.size(); i++)
      f << fields[c].values[i] << " ";
    f << "\n        </DataArray>\n";
  }
  f << "      </CellData>\n";

  f << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
  f.close();
  return true;
}
