#include "fixtures.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/wait.h>

#include "hdf5.h"

namespace fs = std::filesystem;

namespace cardiax_test
{

// ---- meshes ----------------------------------------------------------------

namespace
{

void write_nodes(std::ostream & os, const std::vector<arma::vec3> & pts)
{
  os << "  <nodes size=\"" << pts.size() << "\">\n";
  os << std::setprecision(16);
  for (std::size_t i = 0; i < pts.size(); i++)
    os << "    <node id=\"" << i << "\" x=\"" << pts[i](0) << "\" y=\"" << pts[i](1)
       << "\" z=\"" << pts[i](2) << "\" />\n";
  os << "  </nodes>\n";
}

template <typename Conn>
void write_elements(std::ostream & os, const Conn & elems)
{
  os << "  <elements size=\"" << elems.size() << "\">\n";
  for (std::size_t e = 0; e < elems.size(); e++)
  {
    os << "    <element id=\"" << e << "\"";
    for (std::size_t v = 0; v < elems[e].size(); v++)
      os << " v" << v << "=\"" << elems[e][v] << "\"";
    os << " />\n";
  }
  os << "  </elements>\n";
}

void write_fibers(std::ostream & os, std::size_t nelem, const BoxSpec & b)
{
  const arma::vec3 n = arma::cross(b.fiber, b.sheet);
  os << "  <element_data type=\"fiber_orthotropic\">\n" << std::setprecision(16);
  for (std::size_t e = 0; e < nelem; e++)
  {
    os << "    <element id=\"" << e << "\">"
       << "<fiber>" << b.fiber(0) << "," << b.fiber(1) << "," << b.fiber(2) << "</fiber>"
       << "<sheet>" << b.sheet(0) << "," << b.sheet(1) << "," << b.sheet(2) << "</sheet>"
       << "<normal>" << n(0) << "," << n(1) << "," << n(2) << "</normal>"
       << "</element>\n";
  }
  os << "  </element_data>\n";
}

struct BoundaryFace
{
  int marker;
  std::vector<int> nodes;
};

void write_boundary(std::ostream & os, const std::string & celltype, int dim,
                    const std::vector<BoundaryFace> & faces)
{
  os << "  <boundary celltype=\"" << celltype << "\" dim=\"" << dim << "\">\n";
  for (std::size_t f = 0; f < faces.size(); f++)
  {
    os << "    <element id=\"" << f << "\" marker=\"" << faces[f].marker << "\"";
    for (std::size_t v = 0; v < faces[f].nodes.size(); v++)
      os << " v" << v << "=\"" << faces[f].nodes[v] << "\"";
    os << " />\n";
  }
  os << "  </boundary>\n";
}

std::vector<arma::vec3> box_points(const BoxSpec & b)
{
  std::vector<arma::vec3> pts;
  pts.reserve(box_num_nodes(b));
  for (int k = 0; k <= b.nz; k++)
    for (int j = 0; j <= b.ny; j++)
      for (int i = 0; i <= b.nx; i++)
        pts.push_back(box_point(b, i, j, k));
  return pts;
}

// Outward oriented quads of the box boundary
std::vector<BoundaryFace> box_quads(const BoxSpec & b)
{
  std::vector<BoundaryFace> q;
  auto N = [&](int i, int j, int k) { return box_node(b, i, j, k); };
  for (int k = 0; k < b.nz; k++)
    for (int j = 0; j < b.ny; j++)
    {
      q.push_back({1, {N(0,j,k), N(0,j,k+1), N(0,j+1,k+1), N(0,j+1,k)}});
      q.push_back({2, {N(b.nx,j,k), N(b.nx,j+1,k), N(b.nx,j+1,k+1), N(b.nx,j,k+1)}});
    }
  for (int k = 0; k < b.nz; k++)
    for (int i = 0; i < b.nx; i++)
    {
      q.push_back({3, {N(i,0,k), N(i+1,0,k), N(i+1,0,k+1), N(i,0,k+1)}});
      q.push_back({4, {N(i,b.ny,k), N(i,b.ny,k+1), N(i+1,b.ny,k+1), N(i+1,b.ny,k)}});
    }
  for (int j = 0; j < b.ny; j++)
    for (int i = 0; i < b.nx; i++)
    {
      q.push_back({5, {N(i,j,0), N(i,j+1,0), N(i+1,j+1,0), N(i+1,j,0)}});
      q.push_back({6, {N(i,j,b.nz), N(i+1,j,b.nz), N(i+1,j+1,b.nz), N(i,j+1,b.nz)}});
    }
  return q;
}

} // namespace

void write_tet_box(const std::string & path, const BoxSpec & b)
{
  std::vector<std::vector<int>> tets;
  auto N = [&](int i, int j, int k) { return box_node(b, i, j, k); };
  for (int k = 0; k < b.nz; k++)
    for (int j = 0; j < b.ny; j++)
      for (int i = 0; i < b.nx; i++)
      {
        const int v0 = N(i,j,k),     v1 = N(i+1,j,k),   v2 = N(i+1,j+1,k), v3 = N(i,j+1,k);
        const int v4 = N(i,j,k+1),   v5 = N(i+1,j,k+1), v6 = N(i+1,j+1,k+1), v7 = N(i,j+1,k+1);
        // Kuhn subdivision around the diagonal v0-v6: conforming across cells,
        // all six positively oriented
        tets.push_back({v0, v1, v2, v6});
        tets.push_back({v0, v2, v3, v6});
        tets.push_back({v0, v3, v7, v6});
        tets.push_back({v0, v7, v4, v6});
        tets.push_back({v0, v4, v5, v6});
        tets.push_back({v0, v5, v1, v6});
      }

  // The Kuhn diagonal of each quad face is the one starting at the node with the
  // lowest (i,j,k), which is the first or third node of every quad in box_quads()
  std::vector<BoundaryFace> tris;
  for (const BoundaryFace & q : box_quads(b))
  {
    const std::vector<int> & n = q.nodes;
    const int mn = std::min({n[0], n[1], n[2], n[3]});
    if (n[0] == mn || n[2] == mn)
    {
      tris.push_back({q.marker, {n[0], n[1], n[2]}});
      tris.push_back({q.marker, {n[0], n[2], n[3]}});
    }
    else
    {
      tris.push_back({q.marker, {n[1], n[2], n[3]}});
      tris.push_back({q.marker, {n[1], n[3], n[0]}});
    }
  }

  std::ofstream os(path);
  os << "<?xml version=\"1.0\"?>\n" << b.prefix_xml;
  os << "<mesh celltype=\"tetrahedron\" dim=\"3\">\n";
  write_nodes(os, box_points(b));
  write_elements(os, tets);
  write_fibers(os, tets.size(), b);
  write_boundary(os, "triangle", 2, tris);
  os << b.extra_mesh_xml << "</mesh>\n";
}

void write_hex_box(const std::string & path, const BoxSpec & b)
{
  std::vector<std::vector<int>> hexes;
  auto N = [&](int i, int j, int k) { return box_node(b, i, j, k); };
  for (int k = 0; k < b.nz; k++)
    for (int j = 0; j < b.ny; j++)
      for (int i = 0; i < b.nx; i++)
        hexes.push_back({N(i,j,k), N(i+1,j,k), N(i+1,j+1,k), N(i,j+1,k),
                         N(i,j,k+1), N(i+1,j,k+1), N(i+1,j+1,k+1), N(i,j+1,k+1)});

  std::ofstream os(path);
  os << "<?xml version=\"1.0\"?>\n" << b.prefix_xml;
  os << "<mesh celltype=\"hexahedron\" dim=\"3\">\n";
  write_nodes(os, box_points(b));
  write_elements(os, hexes);
  write_fibers(os, hexes.size(), b);
  write_boundary(os, "quadrilateral", 2, box_quads(b));
  os << b.extra_mesh_xml << "</mesh>\n";
}

namespace
{

std::vector<arma::vec3> square_points(const BoxSpec & b)
{
  std::vector<arma::vec3> pts;
  for (int j = 0; j <= b.ny; j++)
    for (int i = 0; i <= b.nx; i++)
      pts.push_back(arma::vec3({i * b.Lx / b.nx, j * b.Ly / b.ny, 0.0}));
  return pts;
}

std::vector<BoundaryFace> square_edges(const BoxSpec & b)
{
  auto N = [&](int i, int j) { return i + j * (b.nx + 1); };
  std::vector<BoundaryFace> e;
  for (int j = 0; j < b.ny; j++)
  {
    e.push_back({1, {N(0,j+1), N(0,j)}});
    e.push_back({2, {N(b.nx,j), N(b.nx,j+1)}});
  }
  for (int i = 0; i < b.nx; i++)
  {
    e.push_back({3, {N(i,0), N(i+1,0)}});
    e.push_back({4, {N(i+1,b.ny), N(i,b.ny)}});
  }
  return e;
}

} // namespace

void write_tri_square(const std::string & path, const BoxSpec & b)
{
  auto N = [&](int i, int j) { return i + j * (b.nx + 1); };
  std::vector<std::vector<int>> tris;
  for (int j = 0; j < b.ny; j++)
    for (int i = 0; i < b.nx; i++)
    {
      tris.push_back({N(i,j), N(i+1,j), N(i+1,j+1)});
      tris.push_back({N(i,j), N(i+1,j+1), N(i,j+1)});
    }

  std::ofstream os(path);
  os << "<?xml version=\"1.0\"?>\n" << b.prefix_xml;
  os << "<mesh celltype=\"triangle\" dim=\"2\">\n";
  write_nodes(os, square_points(b));
  write_elements(os, tris);
  write_boundary(os, "line", 1, square_edges(b));
  os << b.extra_mesh_xml << "</mesh>\n";
}

void write_quad_square(const std::string & path, const BoxSpec & b)
{
  auto N = [&](int i, int j) { return i + j * (b.nx + 1); };
  std::vector<std::vector<int>> quads;
  for (int j = 0; j < b.ny; j++)
    for (int i = 0; i < b.nx; i++)
      quads.push_back({N(i,j), N(i+1,j), N(i+1,j+1), N(i,j+1)});

  std::ofstream os(path);
  os << "<?xml version=\"1.0\"?>\n" << b.prefix_xml;
  os << "<mesh celltype=\"quadrilateral\" dim=\"2\">\n";
  write_nodes(os, square_points(b));
  write_elements(os, quads);
  write_boundary(os, "line", 1, square_edges(b));
  os << b.extra_mesh_xml << "</mesh>\n";
}

void write_cable(const std::string & path, int n, double L)
{
  std::vector<arma::vec3> pts;
  std::vector<std::vector<int>> segs;
  for (int i = 0; i <= n; i++) pts.push_back(arma::vec3({i * L / n, 0.0, 0.0}));
  for (int i = 0; i < n; i++)  segs.push_back({i, i + 1});

  std::ofstream os(path);
  os << "<?xml version=\"1.0\"?>\n<mesh celltype=\"line\" dim=\"1\">\n";
  write_nodes(os, pts);
  write_elements(os, segs);
  os << "</mesh>\n";
}

std::string fix_nodes_xml(const std::vector<int> & nodes, const std::vector<int> & dirs)
{
  std::ostringstream os;
  for (int n : nodes)
    for (int d : dirs)
      os << "    <node id=\"" << n << "\" direction=\"" << d << "\" value=\"0.0\" />\n";
  return os.str();
}

void write_text(const std::string & path, const std::string & contents)
{
  std::ofstream os(path);
  os << contents;
}

std::string read_text(const std::string & path)
{
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// ---- running the solvers ---------------------------------------------------

AppResult run_app(const std::string & app, const std::string & args, const std::string & dir)
{
  fs::create_directories(dir);
  const std::string logname = "run_" + app + ".log";

  std::ostringstream cmd;
  // A diverged Newton can leave PETSc hanging in its signal handler: never wait forever
  cmd << "cd '" << dir << "' && timeout -k 5 500 '" << app_dir() << "/" << app << "' " << args
      << " > '" << logname << "' 2>&1";

  AppResult r;
  r.dir = dir;
  const int status = std::system(cmd.str().c_str());
  r.exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : 128;
  r.log = read_text(dir + "/" + logname);
  return r;
}

std::string repo_file(const std::string & relative)
{
  return source_dir() + "/" + relative;
}

std::string copy_input(const std::string & relative, const std::string & dir)
{
  fs::create_directories(dir);
  const fs::path src = repo_file(relative);
  fs::copy_file(src, fs::path(dir) / src.filename(), fs::copy_options::overwrite_existing);
  return src.filename().string();
}

double number_after(const std::string & text, const std::string & key)
{
  const std::size_t p = text.rfind(key);
  if (p == std::string::npos) return std::nan("");
  std::istringstream is(text.substr(p + key.size(), 64));
  double v;
  if (is >> v) return v;
  return std::nan("");
}

std::vector<double> numbers_after(const std::string & text, const std::string & key)
{
  std::vector<double> out;
  std::size_t p = text.find(key);
  while (p != std::string::npos)
  {
    std::istringstream is(text.substr(p + key.size(), 64));
    double v;
    if (is >> v) out.push_back(v);
    p = text.find(key, p + key.size());
  }
  return out;
}

arma::mat read_table(const std::string & path)
{
  std::ifstream in(path);
  std::vector<std::vector<double>> rows;
  std::string line;
  while (std::getline(in, line))
  {
    std::istringstream is(line);
    std::vector<double> row;
    double v;
    while (is >> v) row.push_back(v);
    if (!row.empty() && is.eof()) rows.push_back(row);
  }
  if (rows.empty()) return arma::mat();
  arma::mat m(rows.size(), rows[0].size(), arma::fill::zeros);
  for (std::size_t i = 0; i < rows.size(); i++)
    for (std::size_t j = 0; j < rows[i].size() && j < m.n_cols; j++)
      m(i, j) = rows[i][j];
  return m;
}

// ---- HDF5 ------------------------------------------------------------------

namespace
{

herr_t collect_dataset(hid_t, const char * name, const H5O_info_t * info,
                       void * data)
{
  if (info->type == H5O_TYPE_DATASET)
    static_cast<std::vector<std::string> *>(data)->push_back(std::string("/") + name);
  return 0;
}

} // namespace

std::vector<std::string> h5_datasets(const std::string & file)
{
  std::vector<std::string> names;
  hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f < 0) throw std::runtime_error("cannot open HDF5 file " + file);
  H5Ovisit(f, H5_INDEX_NAME, H5_ITER_NATIVE, collect_dataset, &names);
  H5Fclose(f);
  return names;
}

bool h5_has(const std::string & file, const std::string & dataset)
{
  const std::vector<std::string> all = h5_datasets(file);
  const std::string wanted = dataset[0] == '/' ? dataset : "/" + dataset;
  return std::find(all.begin(), all.end(), wanted) != all.end();
}

arma::mat h5_read(const std::string & file, const std::string & dataset)
{
  hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f < 0) throw std::runtime_error("cannot open HDF5 file " + file);
  hid_t d = H5Dopen(f, dataset.c_str(), H5P_DEFAULT);
  if (d < 0)
  {
    H5Fclose(f);
    throw std::runtime_error("no dataset " + dataset + " in " + file);
  }
  hid_t space = H5Dget_space(d);
  const int rank = H5Sget_simple_extent_ndims(space);
  std::vector<hsize_t> dims(std::max(rank, 1), 1);
  if (rank > 0) H5Sget_simple_extent_dims(space, dims.data(), nullptr);

  hsize_t total = 1;
  for (hsize_t n : dims) total *= n;
  std::vector<double> buf(total);
  H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  H5Sclose(space);
  H5Dclose(d);
  H5Fclose(f);

  const arma::uword rows = dims[0];
  const arma::uword cols = total / std::max<hsize_t>(rows, 1);
  arma::mat m(rows, cols);
  for (arma::uword i = 0; i < rows; i++)
    for (arma::uword j = 0; j < cols; j++)
      m(i, j) = buf[i * cols + j];   // HDF5 is row major
  return m;
}

bool h5_equal(const std::string & a, const std::string & b, double rtol,
              std::string * why, const std::vector<std::string> & ignore)
{
  auto fail = [&](const std::string & w) { if (why) *why = w; return false; };

  std::vector<std::string> da = h5_datasets(a), db = h5_datasets(b);
  std::sort(da.begin(), da.end());
  std::sort(db.begin(), db.end());
  if (da != db) return fail("different dataset lists");

  for (const std::string & name : da)
  {
    if (std::find(ignore.begin(), ignore.end(), name) != ignore.end()) continue;
    const arma::mat ma = h5_read(a, name), mb = h5_read(b, name);
    if (ma.n_rows != mb.n_rows || ma.n_cols != mb.n_cols)
      return fail(name + ": different shape");
    for (arma::uword i = 0; i < ma.n_elem; i++)
    {
      if (std::isnan(ma(i)) && std::isnan(mb(i))) continue;
      if (std::isnan(ma(i)) != std::isnan(mb(i))) return fail(name + ": NaN on one side only");
      if (rel_gap(ma(i), mb(i)) > rtol && std::fabs(ma(i) - mb(i)) > 1e-300)
      {
        std::ostringstream os;
        os.precision(16);
        os << name << "[" << i << "]: " << ma(i) << " vs " << mb(i);
        return fail(os.str());
      }
    }
  }
  return true;
}

// ---- regression baseline ---------------------------------------------------

namespace
{

std::string golden_path() { return repo_file("tests/golden/baseline.txt"); }

std::map<std::string, double> read_golden()
{
  std::map<std::string, double> g;
  std::ifstream in(golden_path());
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream is(line);
    std::string k;
    double v;
    if (is >> k >> v) g[k] = v;
  }
  return g;
}

void write_golden(const std::map<std::string, double> & g)
{
  std::ofstream os(golden_path());
  os << "# Regression baseline: values the code produced when the baseline was taken.\n"
     << "# They detect unintended changes and do not prove physical correctness.\n"
     << "# After an INTENDED change: ./run_tests.sh --update-golden, then review the\n"
     << "# diff of this file before keeping it.\n";
  os << std::setprecision(14);
  for (const auto & kv : g) os << kv.first << " " << kv.second << "\n";
}

} // namespace

void check_golden(const char * file, int line, const std::string & key,
                  double value, double rtol)
{
  if (!std::isfinite(value))
  {
    report_failure(file, line, "golden '" + key + "': value is not finite");
    return;
  }

  const char * update = std::getenv("CARDIAX_UPDATE_GOLDEN");
  std::map<std::string, double> g = read_golden();

  if (update && std::string(update) == "1")
  {
    g[key] = value;
    write_golden(g);
    return;
  }

  const auto it = g.find(key);
  if (it == g.end())
  {
    report_failure(file, line, "golden '" + key + "' is missing from tests/golden/baseline.txt "
                               "(got " + show(value) + "); run ./run_tests.sh --update-golden");
    return;
  }
  if (rel_gap(value, it->second) > rtol)
    report_failure(file, line, "golden '" + key + "': got " + show(value) + ", baseline " +
                               show(it->second) + " (rel gap " +
                               show(rel_gap(value, it->second)) + ", tol " + show(rtol) + ")");
}

} // namespace cardiax_test
