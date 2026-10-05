// =============================================================================
//  test_kdtree_uvc.cpp
//
//  Small loop test of the kdtree_uvc.hpp interface: builds the structure ONCE
//  (kdtree_uvc_build) and transfers the field SEVERAL times
//  (kdtree_uvc_transfer), timing the two parts separately.
//
//  It shows what the split is worth: the build pays for the k-NN once and
//  every step of the loop becomes a weighted mean over ready-made lists.
//
//  If the field is a time series, each iteration reads a different step of
//  the HDF5 file (step % n_steps). If it is static, the same vector is
//  reused -- the measured cost is still that of the transfer.
//
//  Usage:
//      ./test_kdtree_uvc Patient_1 Patient_2 --field lat --steps 20
//      ./test_kdtree_uvc P1 P2 --field tecido --type categorical --steps 5
//      ./test_kdtree_uvc P1 P2 --field cell_field/tecido --type categorical
//      ./test_kdtree_uvc P1 P2 --field lat --output target_lat.vtu --nodes 783,5510
// =============================================================================

#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "kdtree_uvc.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// =============================================================================
//  Command line
// =============================================================================

struct Options
{
  string source;
  string target;
  string field;
  string type;
  string output;
  string nodes;         //!< target ids to print, comma separated
  int    steps;
  UVCParameters uvc;

  Options() : source(), target(), field(), type("continuous"), output(),
              nodes(), steps(10), uvc() {}
};

static void help(const char * prog)
{
  cout
    << "usage: " << prog << " <source> <target> --field NAME [options]\n\n"
    << "  source/target   basename or .xmf/.xdmf/.h5 file with nodal UVC\n\n"
    << "options:\n"
    << "  --field NAME    field to transfer (required): just the name\n"
    << "                  (tecido) or the path (vertex_field/tecido,\n"
    << "                  cell_field/tecido) when the name is ambiguous\n"
    << "  --type T        continuous | categorical   (default continuous)\n"
    << "  --steps N       number of kdtree_uvc_transfer calls (default 10)\n"
    << "  --output FILE   writes the result of the last call to a .vtu\n"
    << "  --nodes A,B,C   prints the transferred value at these target ids\n"
    << "  --k N           neighbours per target node (default 12)\n"
    << "  --weight W      gauss | idw            (default gauss)\n"
    << "  --idw-power V   IDW exponent           (default 2)\n"
    << "  --w-ab V --w-tm V --w-rt V   embedding axis weights\n"
    << "  --tv-split V    LV/RV threshold (default: midpoint of the tv range)\n"
    << "  --no-clamp-ab   do not clamp target ab to the source range\n"
    << endl;
}

static bool parse_options(int argc, char ** argv, Options & o)
{
  vector<string> positional;

  for (int i = 1; i < argc; i++)
  {
    const string a = argv[i];

    if (a == "-h" || a == "--help") { help(argv[0]); std::exit(0); }
    else if (a == "--no-clamp-ab") o.uvc.no_clamp_ab = true;
    else if (a.size() > 1 && a[0] == '-')
    {
      if (++i >= argc) { cerr << "[ERROR] " << a << " requires a value." << endl; return false; }
      const string v = argv[i];

      if      (a == "--field")     o.field = v;
      else if (a == "--type")      o.type = v;
      else if (a == "--output")    o.output = v;
      else if (a == "--nodes")     o.nodes = v;
      else if (a == "--steps")     o.steps = atoi(v.c_str());
      else if (a == "--k")         o.uvc.k = atoi(v.c_str());
      else if (a == "--weight")    o.uvc.weight = v;
      else if (a == "--idw-power") o.uvc.idw_power = atof(v.c_str());
      else if (a == "--w-ab")      o.uvc.w_ab = atof(v.c_str());
      else if (a == "--w-tm")      o.uvc.w_tm = atof(v.c_str());
      else if (a == "--w-rt")      o.uvc.w_rt = atof(v.c_str());
      else if (a == "--tv-split")  { o.uvc.tv_split = atof(v.c_str()); o.uvc.has_tv_split = true; }
      else { cerr << "[ERROR] unknown option: " << a << endl; return false; }
    }
    else positional.push_back(a);
  }

  if (positional.size() != 2) { help(argv[0]); return false; }

  o.source = positional[0];
  o.target = positional[1];

  if (o.field.empty())
  {
    cerr << "[ERROR] --field is required." << endl;
    return false;
  }
  if (o.type != "continuous" && o.type != "categorical")
  {
    cerr << "[ERROR] --type must be 'continuous' or 'categorical'." << endl;
    return false;
  }
  if (o.steps < 1) o.steps = 1;

  return true;
}

//! "783,5510,7133" -> {783, 5510, 7133}
static vector<int> parse_ids(const string & s)
{
  vector<int> ids;
  size_t a = 0;
  while (a <= s.size())
  {
    const size_t b = s.find(',', a);
    const string t = s.substr(a, (b == string::npos) ? string::npos : b - a);
    if (!t.empty()) ids.push_back(atoi(t.c_str()));
    if (b == string::npos) break;
    a = b + 1;
  }
  return ids;
}

static double seconds(clock_t a, clock_t b)
{
  return (double) (b - a) / (double) CLOCKS_PER_SEC;
}

// =============================================================================
//  Main program
// =============================================================================

int main(int argc, char ** argv)
{
  Options o;
  if (!parse_options(argc, argv, o)) return 1;

  const bool categorical = (o.type == "categorical");
  cout << std::fixed << std::setprecision(6);
  cout << string(64, '=') << endl;

  // ------------------------------------------------------------ 0) reading
  cout << "Reading meshes..." << endl;

  ReaderHDF5 reader_s, reader_t;
  UVCDataTransfer ms, mt;

  if (!load_uvc_mesh(reader_s, o.source, "source", ms)) return 1;
  if (!load_uvc_mesh(reader_t, o.target, "target", mt)) return 1;

  // ------------------------------------------------------ 1) field (step 0)
  vector<double> source_field;
  bool cell_data = false;
  string path;

  if (!read_field_uvc(reader_s, o.field, source_field, cell_data, path))
  {
    cerr << "[ERROR] field '" << o.field << "' not found in the source." << endl;
    cerr << "        available fields:";
    for (int i = 0; i < reader_s.get_n_fields(); i++)
      cerr << " " << reader_s.get_field(i).path.substr(1)
           << (reader_s.get_field(i).cell_centered ? "(cell)" : "");
    cerr << endl;
    return 1;
  }

  // path is the dataset path ("/vertex_field/tecido"): unique key for the
  // reads in the loop. short_name ("tecido") names the output.
  const int idx = reader_s.find_field(path);
  const int n_steps = (idx >= 0) ? reader_s.get_field(idx).n_steps : 1;
  const string short_name = (idx >= 0) ? reader_s.get_field(idx).name : path;

  cout << "\nField '" << path << "': "
       << (cell_data ? "CellData" : "PointData")
       << ", " << n_steps << " step(s), type " << o.type << endl;

  // ------------------------------------------------------- 2) BUILD (once)
  cout << "\n--- kdtree_uvc_build (1x) ---" << endl;

  KdtreeUVC kd;
  const clock_t t0 = clock();
  if (!kdtree_uvc_build(ms, mt, o.uvc, kd)) return 1;
  const clock_t t1 = clock();

  cout << "  build time: " << seconds(t0, t1) << " s"
       << "   (k_max=" << kd.k_max << ", "
       << kd.copy_to.size() << " copies, "
       << kd.n_unfilled << " unfilled)" << endl;

  // ---------------------------------------------- 3) TRANSFER (in the loop)
  cout << "\n--- kdtree_uvc_transfer (" << o.steps << "x) ---" << endl;

  vector<double> target_field;
  vector<double> raw;
  double t_transfer = 0.0;

  for (int step = 0; step < o.steps; step++)
  {
    // time series: each iteration uses a different step of the file
    if (n_steps > 1)
    {
      const int s = step % n_steps;
      if (!reader_s.read_field_step(path, s, raw)) return 1;
      source_field = cell_data ? cell_to_node_uvc(ms, raw, categorical) : raw;
    }
    else if (step == 0 && cell_data)
    {
      source_field = cell_to_node_uvc(ms, source_field, categorical);
    }

    if ((int) source_field.size() != ms.n_points)
    {
      cerr << "[ERROR] field has " << source_field.size()
           << " values; the source has " << ms.n_points << " nodes." << endl;
      return 1;
    }

    const clock_t a = clock();
    if (!kdtree_uvc_transfer(kd, source_field, categorical, target_field))
      return 1;
    const clock_t b = clock();
    t_transfer += seconds(a, b);

    // --- summary of this iteration ---
    double lo = 1e300, hi = -1e300;
    int n_nan = 0;
    for (int i = 0; i < mt.n_points; i++)
    {
      const double v = target_field[(size_t) i];
      if (!(v == v)) { n_nan++; continue; }
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }

    cout << "  step " << std::setw(3) << step
         << "  output [" << lo << ", " << hi << "]"
         << "  NaN=" << n_nan << endl;
  }

  cout << "\n  total time of the " << o.steps << " transfers: "
       << t_transfer << " s" << endl;
  cout << "  mean per call                  : "
       << t_transfer / (double) o.steps << " s" << endl;
  if (t_transfer > 0.0)
    cout << "  build / mean per call          : "
         << seconds(t0, t1) / (t_transfer / (double) o.steps) << "x" << endl;

  // ------------------------------------------- 4) nodes requested by --nodes
  if (!o.nodes.empty())
  {
    const vector<int> ids = parse_ids(o.nodes);
    cout << "\n--- requested target nodes ---" << endl;
    for (size_t t = 0; t < ids.size(); t++)
    {
      const int i = ids[t];
      if (i < 0 || i >= mt.n_points)
      {
        cout << "  node " << i << ": out of range [0, "
             << mt.n_points - 1 << "]" << endl;
        continue;
      }
      cout << "  node " << i
           << "  " << short_name << "=" << target_field[(size_t) i]
           << "  ab=" << mt.ab[(size_t) i]
           << "  tm=" << mt.tm[(size_t) i]
           << "  tv=" << mt.tv[(size_t) i]
           << "  rt=" << mt.rt[(size_t) i]
           << "  (" << kd.n_neighbors[(size_t) i] << " neighbours)" << endl;
    }
  }

  // ------------------------------------------------------------ 5) output
  if (!o.output.empty())
  {
    vector<OutputField> fields;

    OutputField f_ab; f_ab.name = "ab"; f_ab.values = mt.ab; fields.push_back(f_ab);
    OutputField f_tm; f_tm.name = "tm"; f_tm.values = mt.tm; fields.push_back(f_tm);
    OutputField f_rt; f_rt.name = "rt"; f_rt.values = mt.rt; fields.push_back(f_rt);
    OutputField f_tv; f_tv.name = "tv"; f_tv.values = mt.tv; fields.push_back(f_tv);

    // the field is written with the SAME centering it was read with
    OutputField f_out;
    f_out.name = short_name;
    if (cell_data)
    {
      f_out.cell_data = true;
      f_out.values = node_to_cell_uvc(mt, target_field, categorical);
    }
    else
    {
      f_out.values = target_field;
    }
    fields.push_back(f_out);

    if (!save_vtu_uvc(o.output, mt, fields)) return 1;
    cout << "\nFile saved: " << o.output << endl;
  }

  reader_s.close();
  reader_t.close();

  cout << string(64, '=') << endl;
  cout << "Done" << endl;
  return 0;
}
