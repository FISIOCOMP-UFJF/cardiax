// =============================================================================
//  ReaderHDF5 -- reads the HDF5 files written by WriterHDF5 (HighFive).
//  Format knowledge only; no Mesh. See mesh_hdf5.cpp for building a Mesh.
// =============================================================================

#include "reader_hdf5.hpp"

#include <cctype>
#include <highfive/highfive.hpp>

using std::cerr;
using std::endl;
using std::string;
using std::vector;

namespace
{
  string lower(const string & s)
  {
    string r = s;
    for (size_t i = 0; i < r.size(); i++)
      r[i] = (char) std::tolower((unsigned char) r[i]);
    return r;
  }

  //! Interpret a field dataset shape for n entities (points or elements):
  //!   time series (n_steps x n [x c])  or  static (n [x c]).
  //! Fills n_steps/n_entities/n_comp of f; false if the shape does not fit.
  bool match_shape(const vector<size_t> & dims, int n, FieldInfo & f)
  {
    if (dims.size() >= 2 && dims.size() <= 3 && (int) dims[1] == n)
    {
      f.n_steps    = (int) dims[0];
      f.n_entities = n;
      f.n_comp     = (dims.size() > 2) ? (int) dims[2] : 1;
      return true;
    }
    if (!dims.empty() && dims.size() <= 2 && (int) dims[0] == n)
    {
      f.n_steps    = 1;
      f.n_entities = n;
      f.n_comp     = (dims.size() > 1) ? (int) dims[1] : 1;
      return true;
    }
    return false;
  }

  //! "name", "name.h5", "name.hdf5", "name.xmf", "name.xdmf" -> "name.h5"
  string h5_path(const string & filename)
  {
    const size_t dot = filename.rfind('.');
    const size_t sep = filename.find_last_of("/\\");
    const bool has_ext = (dot != string::npos) &&
                         (sep == string::npos || dot > sep);
    const string ext = has_ext ? lower(filename.substr(dot)) : "";

    if (ext == ".h5" || ext == ".hdf5") return filename;
    if (ext == ".xmf" || ext == ".xdmf") return filename.substr(0, dot) + ".h5";
    return filename + ".h5";
  }
}

// -----------------------------------------------------------------------------

ReaderHDF5::ReaderHDF5()
  : file(), filename(),
    n_node(0), n_points(0), n_elements(0), n_steps(0),
    coords(), connec(), times(), fields()
{
}

ReaderHDF5::~ReaderHDF5()
{
  close();
}

void ReaderHDF5::close()
{
  file.reset();
  filename.clear();
  n_node = n_points = n_elements = n_steps = 0;
  coords.clear();
  connec.clear();
  times.clear();
  fields.clear();
}

// -----------------------------------------------------------------------------

bool ReaderHDF5::open(const string & name)
{
  close();
  const string path = h5_path(name);

  try
  {
    file.reset(new HighFive::File(path, HighFive::File::ReadOnly));
    filename = path;

    // ---------------------------------------------- geometry (np x 3)
    {
      HighFive::DataSet ds = file->getDataSet("/geometry/coordinates");
      const vector<size_t> d = ds.getDimensions();
      if (d.size() != 2 || d[1] < 1 || d[1] > 3)
        throw HighFive::Exception("/geometry/coordinates must be (np x 3)");

      n_points = (int) d[0];
      vector<double> raw(d[0] * d[1]);
      ds.read_raw(raw.data());

      coords.assign((size_t) 3 * n_points, 0.0);       // pads 2D with z = 0
      for (size_t i = 0; i < d[0]; i++)
        for (size_t c = 0; c < d[1]; c++)
          coords[3 * i + c] = raw[d[1] * i + c];
    }

    // ------------------------------------------- topology (ne x nen)
    {
      HighFive::DataSet ds = file->getDataSet("/topology/connectivity");
      const vector<size_t> d = ds.getDimensions();
      if (d.size() != 2)
        throw HighFive::Exception("/topology/connectivity must be (ne x nen)");

      n_elements = (int) d[0];
      n_node     = (int) d[1];
      connec.resize(d[0] * d[1]);
      ds.read_raw(connec.data());
    }

    // ------------------------------------------------------------ time
    // Read as raw values: WriterHDF5 may store it as (n), (n x 1) or
    // (n x 1 x 1); the shape carries no extra information.
    if (file->exist("/time"))
    {
      HighFive::DataSet ds = file->getDataSet("/time");
      times.resize(ds.getElementCount());
      if (!times.empty()) ds.read_raw(times.data());
    }

    // ---------------------------------------------------------- fields
    scan_group("/vertex_field", false);
    scan_group("/cell_field",   true);
  }
  catch (const HighFive::Exception & ex)
  {
    cerr << "ReaderHDF5: failed to open '" << path << "': " << ex.what()
         << endl;
    close();
    return false;
  }

  return true;
}

// -----------------------------------------------------------------------------

void ReaderHDF5::scan_group(const string & group, bool cell)
{
  if (!file->exist(group)) return;

  HighFive::Group g = file->getGroup(group);
  const vector<string> names = g.listObjectNames();

  // Centering is decided by the SHAPE: WriterHDF5 also stores cell data
  // under "/vertex_field". The group name only breaks the tie when
  // n_points == n_elements.
  const int n_first  = cell ? n_elements : n_points;
  const int n_second = cell ? n_points   : n_elements;

  for (size_t k = 0; k < names.size(); k++)
  {
    if (g.getObjectType(names[k]) != HighFive::ObjectType::Dataset) continue;

    const vector<size_t> dims = g.getDataSet(names[k]).getDimensions();

    FieldInfo f;
    f.name = names[k];
    f.path = group + "/" + names[k];

    if (match_shape(dims, n_first, f))
      f.cell_centered = cell;
    else if (match_shape(dims, n_second, f))
    {
      f.cell_centered = !cell;
      cerr << "ReaderHDF5: '" << f.path << "' has one value per "
           << (cell ? "point" : "element") << "; reading it as "
           << (cell ? "PointData" : "CellData") << " (it belongs in "
           << (cell ? "/vertex_field" : "/cell_field") << ")" << endl;
    }
    else
    {
      cerr << "ReaderHDF5: skipping '" << f.path << "': shape matches neither "
           << n_points << " points nor " << n_elements << " elements" << endl;
      continue;
    }

    if (n_steps < f.n_steps) n_steps = f.n_steps;
    fields.push_back(f);
  }
}

// -----------------------------------------------------------------------------

int ReaderHDF5::find_field(const string & name) const
{
  // ------------------------------------------------ by path: "group/name"
  if (name.find('/') != string::npos)
  {
    const string p = (name[0] == '/') ? name : "/" + name;
    for (size_t i = 0; i < fields.size(); i++)
      if (fields[i].path == p) return (int) i;

    const string l = lower(p);
    for (size_t i = 0; i < fields.size(); i++)
      if (lower(fields[i].path) == l) return (int) i;

    return -1;
  }

  // ------------------------------- by name: exact, then case-insensitive
  vector<int> hits;
  for (size_t i = 0; i < fields.size(); i++)
    if (fields[i].name == name) hits.push_back((int) i);

  if (hits.empty())
  {
    const string l = lower(name);
    for (size_t i = 0; i < fields.size(); i++)
      if (lower(fields[i].name) == l) hits.push_back((int) i);
  }

  if (hits.size() == 1) return hits[0];

  if (hits.size() > 1)
  {
    cerr << "ReaderHDF5: field name '" << name << "' is ambiguous; use one of:";
    for (size_t k = 0; k < hits.size(); k++)
      cerr << " " << fields[(size_t) hits[k]].path.substr(1);
    cerr << endl;
  }
  return -1;
}

// -----------------------------------------------------------------------------

bool ReaderHDF5::read_field_step(const string & name, int step,
                                 vector<double> & out) const
{
  if (!is_open())
  {
    cerr << "ReaderHDF5::read_field_step: no file open" << endl;
    return false;
  }

  const int idx = find_field(name);
  if (idx < 0)
  {
    cerr << "ReaderHDF5::read_field_step: field '" << name << "' not found"
         << endl;
    return false;
  }

  const FieldInfo & f = fields[(size_t) idx];
  if (step < 0 || step >= f.n_steps)
  {
    cerr << "ReaderHDF5::read_field_step: step " << step << " out of range [0, "
         << f.n_steps - 1 << "] for '" << f.name << "'" << endl;
    return false;
  }

  out.assign((size_t) f.n_entities * f.n_comp, 0.0);

  try
  {
    HighFive::DataSet ds = file->getDataSet(f.path);
    const vector<size_t> dims = ds.getDimensions();
    const bool series = (dims.size() >= 2 && (int) dims[1] == f.n_entities);

    if (!series)
    {
      ds.read_raw(out.data());                    // static: whole dataset
      return true;
    }

    // hyperslab: one step of (n_steps x n [x c])
    vector<size_t> offset(dims.size(), 0), count(dims);
    offset[0] = (size_t) step;
    count[0]  = 1;
    ds.select(offset, count).read_raw(out.data());
  }
  catch (const HighFive::Exception & ex)
  {
    cerr << "ReaderHDF5::read_field_step: '" << f.path << "' step " << step
         << ": " << ex.what() << endl;
    return false;
  }

  return true;
}

// -----------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & os, const ReaderHDF5 & r)
{
  os << "--- ReaderHDF5: " << r.filename << " ---" << endl;
  os << " points: " << r.n_points << ", elements: " << r.n_elements
     << " (" << r.n_node << " nodes each), steps: " << r.n_steps << endl;
  for (size_t i = 0; i < r.fields.size(); i++)
  {
    const FieldInfo & f = r.fields[i];
    os << "  " << f.name << "  " << (f.cell_centered ? "Cell" : "Node")
       << "  steps=" << f.n_steps << "  comp=" << f.n_comp
       << "  (" << f.path << ")" << endl;
  }
  return os;
}
