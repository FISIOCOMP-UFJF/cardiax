#ifndef READER_HDF5_HPP
#define READER_HDF5_HPP

#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace HighFive { class File; }

/** Description of one field found in an HDF5 file written by WriterHDF5.

    WriterHDF5 stores each field under "/vertex_field" (one value per node)
    or "/cell_field" (one value per element), as a time series:
    (if a dataset's shape contradicts its group, the shape wins and a
    warning is printed when the file is opened)
        scalar : (n_steps x n_entities)
        vector : (n_steps x n_entities x n_comp)
    Static datasets without the step axis are also accepted:
        scalar : (n_entities)
        vector : (n_entities x n_comp)
*/
struct FieldInfo
{
  std::string name;        //!< dataset name ("vm", "ab", ...)
  std::string path;        //!< full HDF5 path ("/vertex_field/ab")
  bool cell_centered;      //!< true = one value per element (CellData)
  int  n_steps;            //!< 1 for static datasets
  int  n_entities;         //!< n_points or n_elements
  int  n_comp;             //!< 1 = scalar, 3 = vector

  FieldInfo()
    : name(), path(), cell_centered(false),
      n_steps(1), n_entities(0), n_comp(1) {}
};

/** Reader for the HDF5 files written by WriterHDF5.

    Knows the FILE FORMAT (where geometry, topology, time and fields are,
    and how their shapes are laid out) and returns everything as plain
    arrays. It does not know what the data is used for: building a Mesh
    (see mesh_hdf5.hpp), transferring UVC fields, post-processing, ...

    Depends only on HighFive/HDF5 -- not on Mesh.

        ReaderHDF5 r;
        if (!r.open("Patient_1_UVC")) return 1;   // or ".h5" / ".xmf"
        std::cout << r;                           // sizes and fields

        const std::vector<double> & xyz = r.get_coordinates();
        std::vector<double> ab;
        r.read_field_step("ab", 0, ab);           // one value per node
*/
class ReaderHDF5
{
public:

  ReaderHDF5();
  ~ReaderHDF5();

  //! Open "<name>.h5". Accepts "name", "name.h5", "name.hdf5" or
  //! "name.xmf" (the matching .h5 is used). Reads geometry, connectivity,
  //! "/time" and the list of fields. False on error (message in cerr).
  bool open(const std::string & filename);

  //! Close the file and clear everything
  void close();

  bool is_open() const { return file.get() != 0; }

  const std::string & get_filename() const { return filename; }

  // ------------------------------------------------------------- mesh data

  int get_n_points()   const { return n_points; }
  int get_n_elements() const { return n_elements; }
  int get_nen()        const { return n_node; }        //!< nodes per element

  //! Flat coordinates: X[3*i + d]
  const std::vector<double> & get_coordinates()  const { return coords; }

  //! Flat connectivity: C[nen*e + j]
  const std::vector<int>    & get_connectivity() const { return connec; }

  // ------------------------------------------------------------------ time

  int get_n_steps() const { return n_steps; }          //!< max over fields

  //! Time values (empty if there is no "/time" dataset)
  const std::vector<double> & get_time() const { return times; }

  // ---------------------------------------------------------------- fields

  int get_n_fields() const { return (int) fields.size(); }
  const FieldInfo & get_field(int i) const { return fields[(size_t) i]; }

  //! Index of a field, or -1. name is either
  //!   a path : "vertex_field/vm", "/cell_field/tecido" (exact dataset)
  //!   a name : "vm" (must be unique across groups; if it exists in both
  //!            /vertex_field and /cell_field, -1 and a message in cerr)
  //! Exact match first, then case-insensitive.
  int find_field(const std::string & name) const;

  //! Read one time step of a field (name or path, as in find_field):
  //! out gets n_entities*n_comp values.
  //! For static fields only step 0 is valid.
  bool read_field_step(const std::string & name, int step,
                       std::vector<double> & out) const;

  //! Summary of the open file (sizes and fields)
  friend std::ostream & operator<<(std::ostream & os, const ReaderHDF5 & r);

private:

  std::unique_ptr<HighFive::File> file;
  std::string filename;

  int n_node, n_points, n_elements, n_steps;

  std::vector<double>    coords;
  std::vector<int>       connec;
  std::vector<double>    times;
  std::vector<FieldInfo> fields;

  //! Register every dataset of a group ("/vertex_field", "/cell_field")
  void scan_group(const std::string & group, bool cell);

  ReaderHDF5(const ReaderHDF5 &) = delete;
  ReaderHDF5 & operator=(const ReaderHDF5 &) = delete;
};

#endif // READER_HDF5_HPP
