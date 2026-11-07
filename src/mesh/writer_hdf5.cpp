#include "writer_hdf5.hpp"
#include <highfive/highfive.hpp>

using HighFive::File;
using HighFive::DataSpace;

namespace {

// strip directory part: "path/to/base" -> "base"
std::string basename_of(const std::string & file)
{
    std::size_t pos = file.find_last_of("/");
    return (pos == std::string::npos) ? file : file.substr(pos + 1);
}

// create a zero-initialised dataset of the given shape
template <typename T>
void create_dataset(File & f, const std::string & path,
                    const std::vector<std::size_t> & shape)
{
    f.createDataSet<T>(path, DataSpace(shape));
}

} // namespace

// ============================================================================
//  Construction / destruction
// ============================================================================

WriterHDF5::WriterHDF5(Mesh * m) : Writer(m)
{
  // do nothing
}

WriterHDF5::~WriterHDF5()
{
  close();
}

// ============================================================================
//  Field registry
// ============================================================================

void WriterHDF5::add_field(const std::string & name, FieldLoc loc,
                           FieldKind kind, int ncomp)
{
    // overwrite if the same name is registered twice
    for (auto & fs : fields_)
        if (fs.name == name) { fs = FieldSpec{name, loc, kind, ncomp}; return; }
    fields_.push_back(FieldSpec{name, loc, kind, ncomp});
}

void WriterHDF5::clear_fields()
{
    fields_.clear();
}

const FieldSpec * WriterHDF5::find_field(const std::string & name) const
{
    for (const auto & fs : fields_)
        if (fs.name == name) return &fs;
    return nullptr;
}

//
// TEMPORARY backward-compatibility set: reproduces exactly the fields the old
// writer hard-coded, so existing solvers keep working unchanged. Remove the
// open()-time call to this once every app registers its own fields.
//
void WriterHDF5::register_default_fields()
{
    // nodal, transient
    add_field("vm",            FieldLoc::Point, FieldKind::Transient);
    add_field("active_stress", FieldLoc::Point, FieldKind::Transient);
    add_field("lambda_f",      FieldLoc::Point, FieldKind::Transient);
    add_field("lambda_rate",   FieldLoc::Point, FieldKind::Transient);
    add_field("displacements", FieldLoc::Point, FieldKind::Transient, 3);

    // cell, transient
    add_field("stress",      FieldLoc::Cell, FieldKind::Transient);
    add_field("strain",      FieldLoc::Cell, FieldKind::Transient);
    add_field("Ta_applied",  FieldLoc::Cell, FieldKind::Transient);
    add_field("ta_scale",    FieldLoc::Cell, FieldKind::Transient);
    add_field("long_strain", FieldLoc::Cell, FieldKind::Transient);
    add_field("circ_strain", FieldLoc::Cell, FieldKind::Transient);
    add_field("rad_strain",  FieldLoc::Cell, FieldKind::Transient);

    // cell, static (mesh-derived; written once instead of replicated per step)
    add_field("fibrosis",   FieldLoc::Cell, FieldKind::Static);
    add_field("aha_marker", FieldLoc::Cell, FieldKind::Static);
}

// legacy adders kept so old call sites still compile
void WriterHDF5::add_ve()
{
    add_field("ve", FieldLoc::Point, FieldKind::Transient);
}

void WriterHDF5::add_fibers()
{
    // kept for compatibility; body was commented out in the original
}

void WriterHDF5::add_scalar_field(std::string & field_name)
{
    add_field(field_name, FieldLoc::Point, FieldKind::Transient);
}

// ============================================================================
//  Open / close
// ============================================================================

void WriterHDF5::open(const std::string & file, int nsteps, double step,
                      bool bido, bool is_restart)
{
    std::string base = basename_of(file);
    if (is_restart) base = base + "_restarted";
    h5name = base + ".h5";

    nsteps_ = static_cast<std::size_t>(nsteps);
    dt_     = step;

    // bidomain adds the extracellular potential field
    if (bido && !find_field("ve"))
        add_field("ve", FieldLoc::Point, FieldKind::Transient);

    // backward compatibility: if the app registered nothing, use the legacy set
    if (fields_.empty())
        register_default_fields();

    write_hdf5(file, nsteps, step);
    write_xdmf(file, nsteps, step, bido);
}

void WriterHDF5::close()
{
    if (xmf.is_open())
        xmf.close();
}

// ============================================================================
//  HDF5 file creation (geometry/topology/time + every registered field)
// ============================================================================

void WriterHDF5::write_hdf5(const std::string & file, int nsteps, double step)
{
    const std::size_t np = mesh->get_n_points();
    const std::size_t ne = mesh->get_n_elements();
    const std::size_t nn = mesh->get_nen();
    const std::size_t ns = static_cast<std::size_t>(nsteps);

    h5name = basename_of(file) + ".h5";
    // keep the _restarted suffix if open() already set it
    // (open() sets h5name before calling us; recompute only the base path here)
    // NOTE: open() is the single source of truth for h5name; this line is a
    // safety net for direct callers of write_hdf5.

    File f(h5name, File::Truncate);
    f.createGroup("/geometry");
    f.createGroup("/topology");
    f.createGroup("/point_data");
    f.createGroup("/cell_data");

    // time array [nsteps, 1]
    {
        std::vector<std::vector<double>> time(ns, std::vector<double>(1));
        for (std::size_t i = 0; i < ns; ++i) time[i][0] = i * step;
        f.createDataSet("/time", time);
    }

    // geometry coordinates [np, 3]
    {
        const std::vector<arma::vec3> & pts = mesh->get_points();
        std::vector<std::vector<double>> coords(np, std::vector<double>(3));
        for (std::size_t i = 0; i < np; ++i) {
            coords[i][0] = pts[i](0);
            coords[i][1] = pts[i](1);
            coords[i][2] = pts[i](2);
        }
        f.createDataSet("/geometry/coordinates", coords);
    }

    // topology connectivity [ne, nn]
    {
        std::vector<std::vector<int>> connec(ne, std::vector<int>(nn));
        for (std::size_t i = 0; i < ne; ++i) {
            std::vector<int> ptnums;
            mesh->get_element_pt_nums(i, ptnums);
            for (std::size_t j = 0; j < nn; ++j) connec[i][j] = ptnums[j];
        }
        f.createDataSet("/topology/connectivity", connec);
    }

    // create every registered field's dataset
    for (const auto & fs : fields_)
    {
        const std::string grp  = (fs.loc == FieldLoc::Point) ? "/point_data/" : "/cell_data/";
        const std::string path = grp + fs.name;
        const std::size_t n    = (fs.loc == FieldLoc::Point) ? np : ne;

        if (fs.kind == FieldKind::Static)
        {
            if (fs.ncomp == 1) create_dataset<double>(f, path, {n});
            else               create_dataset<double>(f, path, {n, (std::size_t)fs.ncomp});
        }
        else // Transient
        {
            if (fs.ncomp == 1) create_dataset<double>(f, path, {ns, n});
            else               create_dataset<double>(f, path, {ns, n, (std::size_t)fs.ncomp});
        }
    }

    // mesh-derived STATIC fields get their values written once here.
    // (Transient fields are written per-step by the solver.)
    if (const FieldSpec * fib = find_field("fibrosis"))
    {
        (void)fib;
        std::vector<double> v(ne);
        for (std::size_t e = 0; e < ne; ++e)
            v[e] = static_cast<double>(mesh->get_element(e).get_index());
        f.getDataSet("/cell_data/fibrosis").write_raw(v.data());
    }
    if (const FieldSpec * aha = find_field("aha_marker"))
    {
        (void)aha;
        std::vector<double> v(ne);
        for (std::size_t e = 0; e < ne; ++e)
            v[e] = static_cast<double>(mesh->get_element(e).get_aha_num());
        f.getDataSet("/cell_data/aha_marker").write_raw(v.data());
    }
}

// ============================================================================
//  Per-step / static writes
// ============================================================================

void WriterHDF5::write_point_field_step(int step, const double * data,
                                        std::string fieldname)
{
    const std::string path = std::string("/point_data/") + fieldname;

    File f(h5name, File::ReadWrite);
    if (!f.exist(path))
    {
        cout << " Warning: nodal field '" << fieldname
             << "' has no dataset; register it with add_field() before open()."
             << endl;
        return;
    }
    f.getDataSet(path)
     .select({(std::size_t)step, 0}, {1, mesh->get_n_points()})
     .write_raw(data);
}

void WriterHDF5::write_cell_field_step(int step, const double * data,
                                       std::string fieldname)
{
    const std::string path = std::string("/cell_data/") + fieldname;

    File f(h5name, File::ReadWrite);
    if (!f.exist(path))
    {
        cout << " Warning: cell field '" << fieldname
             << "' has no dataset; register it with add_field() before open()."
             << endl;
        return;
    }
    f.getDataSet(path)
     .select({(std::size_t)step, 0}, {1, mesh->get_n_elements()})
     .write_raw(data);
}

void WriterHDF5::write_displ_step(int step, const double * displ)
{
    File f(h5name, File::ReadWrite);
    const std::size_t np = mesh->get_n_points();
    f.getDataSet("/point_data/displacements")
     .select({(std::size_t)step, 0, 0}, {1, np, 3})
     .write_raw(displ);
}

void WriterHDF5::write_vm_step(int step, const double * data)
{
    write_point_field_step(step, data, "vm");
}

void WriterHDF5::write_ve_step(int step, const double * data)
{
    write_point_field_step(step, data, "ve");
}

void WriterHDF5::write_static_field(const double * data, std::string fieldname)
{
    const FieldSpec * fs = find_field(fieldname);
    if (!fs) {
        cout << " Warning: static field '" << fieldname
             << "' is not registered." << endl;
        return;
    }
    const std::string grp  = (fs->loc == FieldLoc::Point) ? "/point_data/" : "/cell_data/";
    const std::string path = grp + fieldname;

    File f(h5name, File::ReadWrite);
    if (!f.exist(path)) {
        cout << " Warning: static field '" << fieldname
             << "' has no dataset." << endl;
        return;
    }
    f.getDataSet(path).write_raw(data);
}

// ============================================================================
//  XDMF emission
// ============================================================================

std::string WriterHDF5::xdmf_topology_type() const
{
    const int nn = mesh->get_nen();
    const int nd = mesh->get_n_dim();
    if (nn == 2) return "Polyline";
    if (nn == 3) return "Triangle";
    if (nn == 4 && nd == 2) return "Quadrilateral";
    if (nn == 4 && nd == 3) return "Tetrahedron";
    if (nn == 8) return "Hexahedron";
    return "Hexahedron"; // sensible default
}

//
// One <Attribute> block. Replaces all the hand-copied blocks in the old code.
//  - scalar vs vector  -> AttributeType + hyperslab rank
//  - Point vs Cell     -> Center + count (np/ne) + /point_data vs /cell_data
//  - Static vs Transient -> direct ref vs time-row hyperslab
//
void WriterHDF5::emit_attribute(std::ostream & os, const FieldSpec & fs,
                                int step, int np, int ne, int nsteps) const
{
    const bool        isVec  = (fs.ncomp > 1);
    const char *      atype  = isVec ? "Vector" : "Scalar";
    const char *      center = (fs.loc == FieldLoc::Point) ? "Node" : "Cell";
    const std::size_t n      = (fs.loc == FieldLoc::Point) ? (std::size_t)np : (std::size_t)ne;
    const std::string grp    = (fs.loc == FieldLoc::Point) ? "/point_data/" : "/cell_data/";
    const std::string path   = h5name + ":" + grp + fs.name;

    os << "            <Attribute Name=\"" << fs.name << "\" "
       << "AttributeType=\"" << atype << "\" Center=\"" << center << "\">\n";

    if (fs.kind == FieldKind::Static)
    {
        os << "              <DataItem Format=\"HDF\" NumberType=\"Double\" Precision=\"8\" "
           << "Dimensions=\"" << n;
        if (isVec) os << " 3";
        os << "\">" << path << "</DataItem>\n";
    }
    else // Transient: hyperslab of row `step`
    {
        if (!isVec)
        {
            os << "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"1 " << n << "\" Type=\"HyperSlab\">\n"
               << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
               << "                    " << step << " 0\n"
               << "                    1 1\n"
               << "                    1 " << n << "\n"
               << "                </DataItem>\n"
               << "                <DataItem Dimensions=\"" << nsteps << " " << n << "\" Format=\"HDF\">"
               << path << "</DataItem>\n"
               << "            </DataItem>\n";
        }
        else
        {
            os << "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"1 " << n << " 3\" Type=\"HyperSlab\">\n"
               << "                <DataItem Dimensions=\"3 3\" Format=\"XML\">\n"
               << "                    " << step << " 0 0\n"
               << "                    1 1 1\n"
               << "                    1 " << n << " 3\n"
               << "                </DataItem>\n"
               << "                <DataItem Dimensions=\"" << nsteps << " " << n << " 3\" Format=\"HDF\">"
               << path << "</DataItem>\n"
               << "            </DataItem>\n";
        }
    }

    os << "            </Attribute>\n";
}

void WriterHDF5::write_xdmf(const std::string & file, int nsteps,
                            double /*step*/, bool /*bido*/)
{
    const int np = mesh->get_n_points();
    const int ne = mesh->get_n_elements();
    const int nn = mesh->get_nen();
    const std::string toptype = xdmf_topology_type();

    const std::string xmf_file = basename_of(file) + ".xmf";

    xmf.open(xmf_file.c_str());
    xmf << "<?xml version=\"1.0\" ?>\n"
        << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
        << "<Xdmf Version=\"2.0\">\n"
        << "  <Domain>\n\n";

    // mesh: topology + geometry, referenced by every time grid
    xmf << "    <Topology TopologyType=\"" << toptype << "\"\n"
        << "        NumberOfElements=\"" << ne << "\">\n"
        << "    <DataItem Format=\"HDF\"\n"
        << "        DataType=\"Int\"\n"
        << "        Dimensions=\"" << ne << " " << nn << "\">"
        << h5name << ":/topology/connectivity\n"
        << "    </DataItem>\n"
        << "    </Topology>\n\n";

    xmf << "    <Geometry GeometryType=\"XYZ\">\n"
        << "        <DataItem Dimensions=\"" << np << " 3\"\n"
        << "            NumberType=\"Double\"\n"
        << "            Precision=\"8\"\n"
        << "            Format=\"HDF\">" << h5name << ":/geometry/coordinates\n"
        << "        </DataItem>\n"
        << "    </Geometry>\n\n";

    // temporal collection
    xmf << "    <Grid Name=\"TimeSeries\"\n"
        << "        GridType=\"Collection\"\n"
        << "        CollectionType=\"Temporal\">\n"
        << "        <Time TimeType=\"List\">\n"
        << "            <DataItem Format=\"HDF\"\n"
        << "                NumberType=\"Double\"\n"
        << "                Dimensions=\"" << nsteps << "\">" << h5name << ":/time\n"
        << "            </DataItem>\n"
        << "        </Time>\n";

    // one grid per time step, each referencing the shared mesh
    for (int i = 0; i < nsteps; i++)
    {
        xmf << "        <Grid Name=\"T" << i << "\" GridType=\"Uniform\">\n"
            << "            <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n"
            << "            <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n";

        for (const auto & fs : fields_)
            emit_attribute(xmf, fs, i, np, ne, nsteps);

        xmf << "        </Grid>\n";
    }

    xmf << "    </Grid>\n"
        << "  </Domain>\n"
        << "</Xdmf>\n";

    xmf.close();
}

// ============================================================================
//  Eikonal LAT (static single-field output) -- unchanged behaviour
// ============================================================================

void WriterHDF5::write_eikonal_lat(const std::string & file, const double * lat_data)
{
    std::string base = basename_of(file);
    std::string h5_file  = base + "_lat.h5";
    std::string xmf_file = base + "_lat.xmf";

    const std::size_t np = mesh->get_n_points();
    const std::size_t ne = mesh->get_n_elements();
    const std::size_t nn = mesh->get_nen();

    {
        File f(h5_file, File::Truncate);
        f.createGroup("/geometry");
        f.createGroup("/topology");
        f.createGroup("/point_data");

        std::vector<arma::vec3> pts = mesh->get_points();
        std::vector<std::vector<double>> coords(np, std::vector<double>(3));
        for (std::size_t i = 0; i < np; ++i) {
            coords[i][0] = pts[i](0);
            coords[i][1] = pts[i](1);
            coords[i][2] = pts[i](2);
        }
        f.createDataSet("/geometry/coordinates", coords);

        std::vector<std::vector<int>> connec(ne, std::vector<int>(nn));
        for (std::size_t i = 0; i < ne; ++i) {
            std::vector<int> ptnums;
            mesh->get_element_pt_nums(i, ptnums);
            for (std::size_t j = 0; j < nn; ++j) connec[i][j] = ptnums[j];
        }
        f.createDataSet("/topology/connectivity", connec);

        std::vector<double> lat(lat_data, lat_data + np);
        f.createDataSet("/point_data/lat", lat);
    }

    const std::string toptype = xdmf_topology_type();

    std::ofstream x(xmf_file.c_str());
    x << "<?xml version=\"1.0\" ?>\n"
      << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
      << "<Xdmf Version=\"2.0\">\n"
      << "  <Domain>\n"
      << "    <Grid Name=\"Mesh\" GridType=\"Uniform\">\n"
      << "      <Topology TopologyType=\"" << toptype << "\" NumberOfElements=\"" << ne << "\">\n"
      << "        <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << ne << " " << nn << "\">"
      << h5_file << ":/topology/connectivity</DataItem>\n"
      << "      </Topology>\n"
      << "      <Geometry GeometryType=\"XYZ\">\n"
      << "        <DataItem Format=\"HDF\" NumberType=\"Double\" Precision=\"8\" Dimensions=\"" << np << " 3\">"
      << h5_file << ":/geometry/coordinates</DataItem>\n"
      << "      </Geometry>\n"
      << "      <Attribute Name=\"LAT\" AttributeType=\"Scalar\" Center=\"Node\">\n"
      << "        <DataItem Format=\"HDF\" Dimensions=\"" << np << "\">"
      << h5_file << ":/point_data/lat</DataItem>\n"
      << "      </Attribute>\n"
      << "    </Grid>\n"
      << "  </Domain>\n"
      << "</Xdmf>\n";
    x.close();
}

// ============================================================================
//  Checkpoints (carried over verbatim from the previous implementation)
// ============================================================================

void WriterHDF5::write_checkpoint(int step, double current_time, const double *vm, const double *state_vars, int num_state_vars)
{
    char time_buf[64];
    std::sprintf(time_buf, "%.2f", current_time);

    std::string chk_filename = "checkpoint_t_" + std::string(time_buf) + "ms.h5";

    hid_t file_id = H5Fcreate(chk_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t space_scalar = H5Screate(H5S_SCALAR);

    hid_t attr_step = H5Acreate(file_id, "step", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_step, H5T_NATIVE_INT, &step);
    H5Aclose(attr_step);

    hid_t attr_time = H5Acreate(file_id, "time", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &current_time);
    H5Aclose(attr_time);

    H5Sclose(space_scalar);

    hid_t group_ep = H5Gcreate2(file_id, "/ep", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group_mech = H5Gcreate2(file_id, "/mechanics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t np = mesh->get_n_points();

    hsize_t dims_vm[1] = { np };
    hid_t dataspace_vm = H5Screate_simple(1, dims_vm, NULL);
    hid_t dataset_vm = H5Dcreate(group_ep, "vm", H5T_NATIVE_DOUBLE, dataspace_vm, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_vm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vm);
    H5Dclose(dataset_vm);
    H5Sclose(dataspace_vm);

    hsize_t dims_sv[2] = { np, (hsize_t)num_state_vars };
    hid_t dataspace_sv = H5Screate_simple(2, dims_sv, NULL);
    hid_t dataset_sv = H5Dcreate(group_ep, "state_variables", H5T_NATIVE_DOUBLE, dataspace_sv, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_sv, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, state_vars);
    H5Dclose(dataset_sv);
    H5Sclose(dataspace_sv);

    H5Gclose(group_ep);
    H5Gclose(group_mech);
    H5Fclose(file_id);
}

void WriterHDF5::read_checkpoint_metadata(const std::string &filename, int &step, double &time, int &num_nodes, int &num_vars)
{
    HighFive::File file(filename, HighFive::File::ReadOnly);
    file.getAttribute("step").read(step);
    file.getAttribute("time").read(time);
    auto dims_vm = file.getDataSet("/ep/vm").getSpace().getDimensions();
    num_nodes = (int)dims_vm[0];
    auto dims_sv = file.getDataSet("/ep/state_variables").getSpace().getDimensions();
    num_vars = (int)dims_sv[1];
}

void WriterHDF5::read_checkpoint_data(const std::string &filename, double *vm, double *state_vars)
{
    HighFive::File file(filename, HighFive::File::ReadOnly);
    file.getDataSet("/ep/vm").read_raw(vm);
    file.getDataSet("/ep/state_variables").read_raw(state_vars);
}

void WriterHDF5::write_mech_checkpoint(int step, double current_time, int load_increment, double load_factor,
                                       const double *x_current, const double *fext0, int num_dofs)
{
    std::string final_filename = "checkpoint_step_" + std::to_string(step) + ".h5";
    std::string temp_filename = "checkpoint_mech_tmp.h5";

    hid_t file_id = H5Fcreate(temp_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    {
        hid_t space_scalar = H5Screate(H5S_SCALAR);
        hid_t attr_step = H5Acreate(file_id, "step", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_step, H5T_NATIVE_INT, &step);
        H5Aclose(attr_step);
        hid_t attr_time = H5Acreate(file_id, "time", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &current_time);
        H5Aclose(attr_time);
        H5Sclose(space_scalar);
    }

    hid_t group_mech;
    if (H5Lexists(file_id, "/mechanics", H5P_DEFAULT) > 0)
        group_mech = H5Gopen2(file_id, "/mechanics", H5P_DEFAULT);
    else
        group_mech = H5Gcreate2(file_id, "/mechanics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t space_scalar = H5Screate(H5S_SCALAR);

    if (H5Aexists(group_mech, "load_increment")) H5Adelete(group_mech, "load_increment");
    hid_t attr_load = H5Acreate(group_mech, "load_increment", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_load, H5T_NATIVE_INT, &load_increment);
    H5Aclose(attr_load);

    if (H5Aexists(group_mech, "load_factor")) H5Adelete(group_mech, "load_factor");
    hid_t attr_lf = H5Acreate(group_mech, "load_factor", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_lf, H5T_NATIVE_DOUBLE, &load_factor);
    H5Aclose(attr_lf);

    H5Sclose(space_scalar);

    auto save_array = [&](const char* name, hid_t space, const double* data) {
        hid_t dset;
        if (H5Lexists(group_mech, name, H5P_DEFAULT) > 0)
            dset = H5Dopen2(group_mech, name, H5P_DEFAULT);
        else
            dset = H5Dcreate(group_mech, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(dset);
    };

    hsize_t dims[1] = { (hsize_t)num_dofs };
    hid_t space_nodal = H5Screate_simple(1, dims, NULL);
    save_array("x_current", space_nodal, x_current);
    save_array("fext0", space_nodal, fext0);
    H5Sclose(space_nodal);
    H5Gclose(group_mech);
    H5Fclose(file_id);

    std::rename(temp_filename.c_str(), final_filename.c_str());
}

void WriterHDF5::read_mech_checkpoint_metadata(const std::string &filename,
                                               int &step, double &time, int &load_increment,
                                               double &load_factor, int &num_dofs)
{
    HighFive::File file(filename, HighFive::File::ReadOnly);
    file.getAttribute("step").read(step);
    file.getAttribute("time").read(time);
    auto group_mech = file.getGroup("/mechanics");
    group_mech.getAttribute("load_increment").read(load_increment);
    group_mech.getAttribute("load_factor").read(load_factor);
    auto dims = group_mech.getDataSet("x_current").getSpace().getDimensions();
    num_dofs = (int)dims[0];
}

void WriterHDF5::read_mech_checkpoint_data(const std::string &filename,
                                           double *x_current, double *fext0)
{
    HighFive::File file(filename, HighFive::File::ReadOnly);
    auto group_mech = file.getGroup("/mechanics");
    group_mech.getDataSet("x_current").read_raw(x_current);
    group_mech.getDataSet("fext0").read_raw(fext0);
}

void WriterHDF5::write_coupled_checkpoint(int step, double current_time,
                                          const double *vm, const double *state_vars, int num_state_vars,
                                          int load_increment, double load_factor,
                                          const double *x_current, const double *fext0, int num_dofs)
{
    char time_buf[64];
    std::sprintf(time_buf, "%.2f", current_time);

    std::string final_filename = "checkpoint_t_" + std::string(time_buf) + "ms.h5";
    std::string temp_filename = "checkpoint_coupled_tmp.h5";

    hid_t file_id = H5Fcreate(temp_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t space_scalar = H5Screate(H5S_SCALAR);

    hid_t attr_step = H5Acreate(file_id, "step", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_step, H5T_NATIVE_INT, &step);
    H5Aclose(attr_step);

    hid_t attr_time = H5Acreate(file_id, "time", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &current_time);
    H5Aclose(attr_time);

    hid_t group_ep = H5Gcreate2(file_id, "/ep", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t np = mesh->get_n_points();

    hsize_t dims_vm[1] = { np };
    hid_t space_vm = H5Screate_simple(1, dims_vm, NULL);
    hid_t dataset_vm = H5Dcreate(group_ep, "vm", H5T_NATIVE_DOUBLE, space_vm, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_vm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vm);
    H5Dclose(dataset_vm);
    H5Sclose(space_vm);

    hsize_t dims_sv[2] = { np, (hsize_t)num_state_vars };
    hid_t space_sv = H5Screate_simple(2, dims_sv, NULL);
    hid_t dataset_sv = H5Dcreate(group_ep, "state_variables", H5T_NATIVE_DOUBLE, space_sv, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_sv, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, state_vars);
    H5Dclose(dataset_sv);
    H5Sclose(space_sv);

    H5Gclose(group_ep);

    hid_t group_mech = H5Gcreate2(file_id, "/mechanics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t attr_load = H5Acreate(group_mech, "load_increment", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_load, H5T_NATIVE_INT, &load_increment);
    H5Aclose(attr_load);

    hid_t attr_lf = H5Acreate(group_mech, "load_factor", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_lf, H5T_NATIVE_DOUBLE, &load_factor);
    H5Aclose(attr_lf);

    hsize_t dims_mech[1] = { (hsize_t)num_dofs };
    hid_t space_mech = H5Screate_simple(1, dims_mech, NULL);

    hid_t dset_x = H5Dcreate(group_mech, "x_current", H5T_NATIVE_DOUBLE, space_mech, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x_current);
    H5Dclose(dset_x);

    hid_t dset_fext = H5Dcreate(group_mech, "fext0", H5T_NATIVE_DOUBLE, space_mech, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_fext, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fext0);
    H5Dclose(dset_fext);

    H5Sclose(space_mech);
    H5Gclose(group_mech);

    H5Sclose(space_scalar);
    H5Fclose(file_id);

    std::rename(temp_filename.c_str(), final_filename.c_str());
}



// #include "writer_hdf5.hpp"

// #include <highfive/highfive.hpp>

// using HighFive::File;
// using HighFive::DataSpace;

// namespace {

// // create a zero-initialised dataset of the given shape
// template <typename T>
// void create_dataset(File & f, const std::string & path,
//                     const std::vector<std::size_t> & shape)
// {
//     f.createDataSet<T>(path, DataSpace(shape));
// }

// // write one time slice (row `step`) into an existing 2-D dataset [nsteps, ncols]
// void write_row(const std::string & h5name, const std::string & path,
//                std::size_t step, const double * data, std::size_t ncols)
// {
//     File f(h5name, File::ReadWrite);
//     auto dset = f.getDataSet(path);
//     dset.select({step, 0}, {1, ncols}).write_raw(data);
// }

// } // namespace

// WriterHDF5::WriterHDF5(Mesh * m) : Writer(m)
// {
//   // do nothing
// }

// WriterHDF5::~WriterHDF5()
// {
//   close();
// }

// void WriterHDF5::open(const std::string & file, int nsteps, double step, bool bido, bool is_restart)
// {
//     std::size_t pos  = file.find_last_of("/");
//     std::string base = file.c_str();
//     if (pos != std::string::npos) {
//         base = file.substr(pos+1);
//     }
//     if (is_restart) base = base + "_restarted";
//     h5name = base + ".h5";

//     write_hdf5(file, nsteps, step);
//     write_xdmf(file, nsteps, step, bido);
// }

// void WriterHDF5::close()
// {
//   // check if need to do something else
// }

// void WriterHDF5::write_hdf5(const std::string & file, int nsteps, double step)
// {
//     const std::size_t np = mesh->get_n_points();
//     const std::size_t ne = mesh->get_n_elements();
//     const std::size_t nn = mesh->get_nen();
//     const std::size_t ns = static_cast<std::size_t>(nsteps);

//     std::size_t pos  = file.find_last_of("/");
//     std::string base = file.c_str();
//     if (pos != std::string::npos)
//       base = file.substr(pos+1);
//     h5name = base + ".h5";

//     // create (truncate) the file and the groups
//     File f(h5name, File::Truncate);
//     f.createGroup("/geometry");
//     f.createGroup("/topology");
//     f.createGroup("/point_data");
//     f.createGroup("/cell_data");

//     // time array [nsteps, 1]
//     {
//         std::vector<std::vector<double>> time(ns, std::vector<double>(1));
//         for (std::size_t i = 0; i < ns; ++i) time[i][0] = i * step;
//         f.createDataSet("/time", time);
//     }

//     // geometry coordinates [np, 3]
//     {
//         const std::vector<arma::vec3> & pts = mesh->get_points();
//         std::vector<std::vector<double>> coords(np, std::vector<double>(3));
//         for (std::size_t i = 0; i < np; ++i) {
//             coords[i][0] = pts[i](0);
//             coords[i][1] = pts[i](1);
//             coords[i][2] = pts[i](2);
//         }
//         f.createDataSet("/geometry/coordinates", coords);
//     }

//     // topology connectivity [ne, nn]
//     // TODO: mesh is fixed for one type of element only -> improve this
//     {
//         std::vector<std::vector<int>> connec(ne, std::vector<int>(nn));
//         for (std::size_t i = 0; i < ne; ++i) {
//             std::vector<int> ptnums;
//             mesh->get_element_pt_nums(i, ptnums);
//             for (std::size_t j = 0; j < nn; ++j) connec[i][j] = ptnums[j];
//         }
//         f.createDataSet("/topology/connectivity", connec);
//     }

//     // nodal fields [nsteps, np]
//     create_dataset<double>(f, "/point_data/vm",            {ns, np});
//     create_dataset<double>(f, "/point_data/active_stress", {ns, np});
//     create_dataset<double>(f, "/point_data/lambda_f",      {ns, np});
//     create_dataset<double>(f, "/point_data/lambda_rate",   {ns, np});

//     // cell fields [nsteps, ne]
//     create_dataset<double>(f, "/cell_data/stress",      {ns, ne});
//     create_dataset<double>(f, "/cell_data/strain",      {ns, ne});
//     create_dataset<double>(f, "/cell_data/Ta_applied",  {ns, ne});
//     create_dataset<double>(f, "/cell_data/ta_scale",    {ns, ne});
//     create_dataset<double>(f, "/cell_data/long_strain", {ns, ne});
//     create_dataset<double>(f, "/cell_data/circ_strain", {ns, ne});
//     create_dataset<double>(f, "/cell_data/rad_strain",  {ns, ne});

//     // fibrosis + aha_marker: per-element values replicated across all steps,
//     // stored as doubles to match the existing on-disk layout
//     {
//         std::vector<std::vector<double>> fib(ns, std::vector<double>(ne));
//         for (std::size_t i = 0; i < ns; ++i)
//             for (std::size_t e = 0; e < ne; ++e)
//                 fib[i][e] = static_cast<double>(mesh->get_element(e).get_index());
//         f.createDataSet("/cell_data/fibrosis", fib);

//         std::vector<std::vector<double>> aha(ns, std::vector<double>(ne));
//         for (std::size_t i = 0; i < ns; ++i)
//             for (std::size_t e = 0; e < ne; ++e)
//                 aha[i][e] = static_cast<double>(mesh->get_element(e).get_aha_num());
//         f.createDataSet("/cell_data/aha_marker", aha);
//     }

//     // displacement vector field [nsteps, np, 3]
//     create_dataset<double>(f, "/point_data/displacements", {ns, np, 3});
// }

// void WriterHDF5::write_cell_field_step(int step, const double *data, string fieldname)
// {
//     write_row(h5name, std::string("/cell_data/") + fieldname,
//               static_cast<std::size_t>(step), data, mesh->get_n_elements());
// }

// void WriterHDF5::write_eikonal_lat(const std::string & file, const double *lat_data)
// {
//     std::string base = file;
//     std::size_t pos = file.find_last_of("/");
//     if (pos != std::string::npos) base = file.substr(pos + 1);

//     std::string h5_file  = base + "_lat.h5";
//     std::string xmf_file = base + "_lat.xmf";

//     const std::size_t np = mesh->get_n_points();
//     const std::size_t ne = mesh->get_n_elements();
//     const std::size_t nn = mesh->get_nen();

//     // hdf5 part
//     {
//         File f(h5_file, File::Truncate);
//         f.createGroup("/geometry");
//         f.createGroup("/topology");
//         f.createGroup("/point_data");

//         std::vector<arma::vec3> pts = mesh->get_points();
//         std::vector<std::vector<double>> coords(np, std::vector<double>(3));
//         for (std::size_t i = 0; i < np; ++i) {
//             coords[i][0] = pts[i](0);
//             coords[i][1] = pts[i](1);
//             coords[i][2] = pts[i](2);
//         }
//         f.createDataSet("/geometry/coordinates", coords);

//         std::vector<std::vector<int>> connec(ne, std::vector<int>(nn));
//         for (std::size_t i = 0; i < ne; ++i) {
//             std::vector<int> ptnums;
//             mesh->get_element_pt_nums(i, ptnums);
//             for (std::size_t j = 0; j < nn; ++j) connec[i][j] = ptnums[j];
//         }
//         f.createDataSet("/topology/connectivity", connec);

//         // lat [np]
//         std::vector<double> lat(lat_data, lat_data + np);
//         f.createDataSet("/point_data/lat", lat);
//     }

//     // xdmf part
//     int nd = mesh->get_n_dim();
//     std::string toptype;
//     if (nn == 2) toptype = "Polyline";
//     else if (nn == 3) toptype = "Triangle";
//     else if (nn == 4 && nd == 2) toptype = "Quadrilateral";
//     else if (nn == 4 && nd == 3) toptype = "Tetrahedron";
//     else if (nn == 8) toptype = "Hexahedron";

//     std::ofstream xmf(xmf_file.c_str());
//     xmf << "<?xml version=\"1.0\" ?>\n"
//         << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
//         << "<Xdmf Version=\"2.0\">\n"
//         << "  <Domain>\n"
//         << "    <Grid Name=\"Mesh\" GridType=\"Uniform\">\n"
//         << "      <Topology TopologyType=\"" << toptype << "\" NumberOfElements=\"" << ne << "\">\n"
//         << "        <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << ne << " " << nn << "\">"
//         << h5_file << ":/topology/connectivity</DataItem>\n"
//         << "      </Topology>\n"
//         << "      <Geometry GeometryType=\"XYZ\">\n"
//         << "        <DataItem Format=\"HDF\" NumberType=\"Double\" Precision=\"8\" Dimensions=\"" << np << " 3\">"
//         << h5_file << ":/geometry/coordinates</DataItem>\n"
//         << "      </Geometry>\n"
//         << "      <Attribute Name=\"LAT\" AttributeType=\"Scalar\" Center=\"Node\">\n"
//         << "        <DataItem Format=\"HDF\" Dimensions=\"" << np << "\">"
//         << h5_file << ":/point_data/lat</DataItem>\n"
//         << "      </Attribute>\n"
//         << "    </Grid>\n"
//         << "  </Domain>\n"
//         << "</Xdmf>\n";
//     xmf.close();
// }

// void WriterHDF5::write_point_field_step(int step, const double *data,
//                                         string fieldname)
// {
//     const std::string path = std::string("/point_data/") + fieldname;

//     File f(h5name, File::ReadWrite);

//     // a field that was never created cannot be written into; warn once and skip
//     if (!f.exist(path))
//     {
//         cout << " Warning: nodal field '" << fieldname
//              << "' has no dataset in the HDF5 file; not written."
//              << " Add it to WriterHDF5::write_hdf5() and"
//              << " WriterHDF5::write_xdmf()." << endl;
//         return;
//     }

//     auto dset = f.getDataSet(path);
//     dset.select({static_cast<std::size_t>(step), 0},
//                 {1, mesh->get_n_points()}).write_raw(data);
// }

// void WriterHDF5::write_vm_step(int step, const double *data)
// {
//     write_row(h5name, "/point_data/vm",
//               static_cast<std::size_t>(step), data, mesh->get_n_points());
// }

// void WriterHDF5::write_displ_step(int step, const double *displ)
// {
//     // displacement [nsteps, np, 3]; write the [1, np, 3] slice for `step`
//     File f(h5name, File::ReadWrite);
//     auto dset = f.getDataSet("/point_data/displacements");
//     const std::size_t np = mesh->get_n_points();
//     dset.select({static_cast<std::size_t>(step), 0, 0}, {1, np, 3}).write_raw(displ);
// }

// void WriterHDF5::add_ve()
// {
//     // create /point_data/ve with the same [nsteps, nnodes] shape as vm
//     File f(h5name, File::ReadWrite);
//     auto vm = f.getDataSet("/point_data/vm");
//     auto dims = vm.getSpace().getDimensions();
//     f.createDataSet<double>("/point_data/ve", DataSpace(dims));
// }

// void WriterHDF5::add_scalar_field(std::string & field_name)
// {
//     // size the new field from an existing field's dimensions; as in the
//     // original, the created dataset is /point_data/ve
//     File f(h5name, File::ReadWrite);
//     auto src = f.getDataSet(std::string("/point_data/") + field_name);
//     auto dims = src.getSpace().getDimensions();
//     f.createDataSet<double>("/point_data/ve", DataSpace(dims));
// }

// void WriterHDF5::add_fibers()
// {
//     // unchanged: body was commented out in the original
// }

// void WriterHDF5::write_ve_step(int step, const double *data)
// {
//     write_row(h5name, "/point_data/ve",
//               static_cast<std::size_t>(step), data, mesh->get_n_points());
// }

// void WriterHDF5::write_xdmf(const std::string & file, int nsteps, 
//                             double step, bool bido)
// {
//     int np = mesh->get_n_points();
//     int ne = mesh->get_n_elements();
//     int nn = mesh->get_nen(); 
//     int nd = mesh->get_n_dim();

//     std::string toptype;
//     if (nn == 2) toptype = "Polyline";
//     else if (nn == 3) toptype = "Triangle";
//     else if (nn == 4 && nd == 2) toptype = "Quadrilateral";
//     else if (nn == 4 && nd == 3) toptype = "Tetrahedron";
//     else if (nn == 8) toptype = "Hexahedron";

//     std::size_t pos  = file.find_last_of("/");
//     std::string base = file.c_str();
//     if (pos != std::string::npos)
//       base = file.substr(pos+1);
//     std::string xmf_file(base + ".xmf");
    
//     if(bido)
//     {
//         // do nothing for now...
//     }
        
//     //
//     // open the file and write the XML description of the mesh
//     //
//     xmf.open(xmf_file.c_str());
//     xmf << "<?xml version=\"1.0\" ?>\n"
//         << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
//         << "<Xdmf Version=\"2.0\">\n"
//         << "  <Domain>\n\n";
            
//     xmf << "    <Topology TopologyType=\"" << toptype << "\"\n"
//         << "        NumberOfElements=\"" << ne << "\">\n"
//         << "    <DataItem Format=\"HDF\" \n"
//         << "        DataType=\"Int\"  \n"
//         << "        Dimensions=\"" << ne << " " << nn << "\">" 
//         << h5name << ":/topology/connectivity\n"
//         << "    </DataItem>\n"
//         << "    </Topology>\n\n";        
    
//     xmf << "    <Geometry GeometryType=\"XYZ\">\n"
//         << "        <DataItem Dimensions=\"" << np << " " << 3 << "\"\n"
//         << "            NumberType=\"Double\" \n"
//         << "            Precision=\"8\" \n"
//         << "            Format=\"HDF\">" << h5name << ":/geometry/coordinates\n"
//         << "        </DataItem>\n"
//         << "    </Geometry>\n\n";        
    
//     //
//     // write time information
//     //    
//     xmf << "    <Grid Name=\"TimeSeries\" \n"
//         << "        GridType=\"Collection\"\n"
//         << "        CollectionType=\"Temporal\">\n"           
//         << "        <Time TimeType=\"List\">\n" // for NON-UNIFORM time steps
//         << "            <DataItem  Format=\"HDF\" \n" 
//         << "                NumberType=\"Double\" \n"
//         << "                Dimensions=\"" << nsteps << "\">" << h5name << ":/time\n"
//         << "            </DataItem>\n"
//         << "        </Time>\n";
    
//     for(int i=0; i<nsteps; i++)
//     {   
//         xmf << "        <Grid Name=\"T" << i << "\" GridType=\"Uniform\">\n"
//             << "            <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n"
//             << "            <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n"; 
        
//         //
//         // potential
//         //     
//         xmf << "            <Attribute Name=\"vm\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Node\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << np << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << np <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Points\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/point_data/vm\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";

//         //
//         // active tension (nodal)
//         //
//         xmf << "            <Attribute Name=\"active_stress\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Node\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << np << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << np <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Points\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/point_data/active_stress\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";

//         //
//         // fibre stretch lambda_f (nodal)
//         //
//         xmf << "            <Attribute Name=\"lambda_f\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Node\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << np << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << np <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Points\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/point_data/lambda_f\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";

//         //
//         // fibre stretch rate d(lambda_f)/dt (nodal)
//         //
//         xmf << "            <Attribute Name=\"lambda_rate\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Node\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << np << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << np <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Points\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/point_data/lambda_rate\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";

//         //teste cell field
//         //Ta_applied -- tensao ativa efetivamente aplicada pela mecanica
//         xmf << "            <Attribute Name=\"Ta_applied\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Cell\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << ne << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << ne <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Cells\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/cell_data/Ta_applied\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//         //ta_scale -- mapa do multiplicador por elemento
//         xmf << "            <Attribute Name=\"ta_scale\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Cell\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << ne << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << ne <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Cells\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/cell_data/ta_scale\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//         //stress
//         xmf << "            <Attribute Name=\"stress\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Cell\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << ne << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << ne <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Cells\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/cell_data/stress\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//         //strain
//         xmf << "            <Attribute Name=\"strain\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Cell\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << ne << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << ne <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Cells\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/cell_data/strain\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//         //long_strain
//         xmf << "            <Attribute Name=\"long_strain\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Cell\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << ne << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << ne <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Cells\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/cell_data/long_strain\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//         //circ_strain
//         xmf << "            <Attribute Name=\"circ_strain\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Cell\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << ne << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << ne <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Cells\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/cell_data/circ_strain\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//         //rad_strain
//         xmf << "            <Attribute Name=\"rad_strain\" \n"
//             << "                AttributeType=\"Scalar\" \n"
//             << "                Center=\"Cell\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << ne << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//             << "                    " << i << " 0 \n"
//             << "                    1 1 \n"
//             << "                    1 " << ne <<"\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Cells\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/cell_data/rad_strain\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//   	//fibrosis
//       xmf << "            <Attribute Name=\"fibrosis\" \n"
//           << "                AttributeType=\"Scalar\" \n"
//           << "                Center=\"Cell\">\n"
//           << "            <DataItem ItemType=\"HyperSlab\" \n"
//           << "                Dimensions=\"1 " << ne << "\" \n"
//           << "                Type=\"HyperSlab\">\n"
//           << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//           << "                    " << i << " 0 \n"
//           << "                    1 1 \n"
//           << "                    1 " << ne <<"\n"
//           << "                </DataItem>\n"
//           << "                <DataItem Name=\"Cells\" \n"
//           << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//           << "                    Format=\"HDF\">" << h5name << ":/cell_data/fibrosis\n"
//           << "                </DataItem>\n"
//           << "            </DataItem>\n"
//           << "            </Attribute>\n";
//     //aha_segment
//       xmf << "            <Attribute Name=\"aha_marker\" \n"
//           << "                AttributeType=\"Scalar\" \n"
//           << "                Center=\"Cell\">\n"
//           << "            <DataItem ItemType=\"HyperSlab\" \n"
//           << "                Dimensions=\"1 " << ne << "\" \n"
//           << "                Type=\"HyperSlab\">\n"
//           << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//           << "                    " << i << " 0 \n"
//           << "                    1 1 \n"
//           << "                    1 " << ne <<"\n"
//           << "                </DataItem>\n"
//           << "                <DataItem Name=\"Cells\" \n"
//           << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
//           << "                    Format=\"HDF\">" << h5name << ":/cell_data/aha_marker\n"
//           << "                </DataItem>\n"
//           << "            </DataItem>\n"
//           << "            </Attribute>\n";
//           //fim teste
        
//         if(bido)
//         {
//             //
//             // extracellular potential
//             //     
//             xmf << "            <Attribute Name=\"ve\" \n"
//                 << "                AttributeType=\"Scalar\" \n"
//                 << "                Center=\"Node\">\n"
//                 << "            <DataItem ItemType=\"HyperSlab\" \n"
//                 << "                Dimensions=\"1 " << np << "\" \n"
//                 << "                Type=\"HyperSlab\">\n"
//                 << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
//                 << "                    " << i << " 0 \n"
//                 << "                    1 1 \n"
//                 << "                    1 " << np <<"\n"
//                 << "                </DataItem>\n"
//                 << "                <DataItem Name=\"Points\" \n"
//                 << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
//                 << "                    Format=\"HDF\">" << h5name << ":/point_data/ve\n"
//                 << "                </DataItem>\n"
//                 << "            </DataItem>\n"
//                 << "            </Attribute>\n";          
//         }
        
//         //
//         // displacements
//         //
//         xmf << "            <Attribute Name=\"displacement\" \n"
//             << "                AttributeType=\"Vector\" \n"
//             << "                Center=\"Node\">\n"
//             << "            <DataItem ItemType=\"HyperSlab\" \n"
//             << "                Dimensions=\"1 " << np  << " " << 3 << "\" \n"
//             << "                Type=\"HyperSlab\">\n"
//             << "                <DataItem Dimensions=\"3 3\" Format=\"XML\">\n"
//             << "                    " << i << " 0 0 \n"
//             << "                    1 1 1 \n"
//             << "                    1 " << np << " 3" << "\n"
//             << "                </DataItem>\n"
//             << "                <DataItem Name=\"Points\" \n"
//             << "                    Dimensions=\"" << nsteps << " " << np << " 3" << "\" \n"
//             << "                    Format=\"HDF\">" << h5name << ":/point_data/displacements\n"
//             << "                </DataItem>\n"
//             << "            </DataItem>\n"
//             << "            </Attribute>\n";
//         xmf << "        </Grid>\n";
//     }
    
//     // Used to write only the mesh
//     //fprintf(xmf, "    <Grid Name=\"Mesh\" GridType=\"Uniform\">\n");
//     //fprintf(xmf, "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n");
//     //fprintf(xmf, "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");
//     //fprintf(xmf, "    </Grid>\n\n");     
    
//     //
//     // close file
//     //
//     xmf << "    </Grid>\n"
//         << "  </Domain>\n"
//         << "</Xdmf>\n";
    
//     xmf.close();
// }

// void WriterHDF5::write_checkpoint(int step, double current_time, const double *vm, const double *state_vars, int num_state_vars)
// {
//     char time_buf[64];
//     std::sprintf(time_buf, "%.2f", current_time); 
    
//     std::string chk_filename = "checkpoint_t_" + std::string(time_buf) + "ms.h5";

//     hid_t file_id = H5Fcreate(chk_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
//     hid_t space_scalar = H5Screate(H5S_SCALAR);
    
//     hid_t attr_step = H5Acreate(file_id, "step", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_step, H5T_NATIVE_INT, &step);
//     H5Aclose(attr_step);

//     hid_t attr_time = H5Acreate(file_id, "time", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &current_time);
//     H5Aclose(attr_time);
    
//     H5Sclose(space_scalar);

//     hid_t group_ep = H5Gcreate2(file_id, "/ep", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     hid_t group_mech = H5Gcreate2(file_id, "/mechanics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

//     hsize_t np = mesh->get_n_points();

//     // saving vm
//     hsize_t dims_vm[1] = { np };
//     hid_t dataspace_vm = H5Screate_simple(1, dims_vm, NULL);
//     hid_t dataset_vm = H5Dcreate(group_ep, "vm", H5T_NATIVE_DOUBLE, dataspace_vm, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
//     H5Dwrite(dataset_vm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vm);
    
//     H5Dclose(dataset_vm);
//     H5Sclose(dataspace_vm);

//     // saving state_variables
//     hsize_t dims_sv[2] = { np, (hsize_t)num_state_vars };
//     hid_t dataspace_sv = H5Screate_simple(2, dims_sv, NULL);
//     hid_t dataset_sv = H5Dcreate(group_ep, "state_variables", H5T_NATIVE_DOUBLE, dataspace_sv, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
//     H5Dwrite(dataset_sv, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, state_vars);
    
//     H5Dclose(dataset_sv);
//     H5Sclose(dataspace_sv);

//     H5Gclose(group_ep);
//     H5Gclose(group_mech);
//     H5Fclose(file_id);
// }

// void WriterHDF5::read_checkpoint_metadata(const std::string &filename, int &step, double &time, int &num_nodes, int &num_vars)
// {
//     HighFive::File file(filename, HighFive::File::ReadOnly);

//     file.getAttribute("step").read(step);
//     file.getAttribute("time").read(time);

//     auto dims_vm = file.getDataSet("/ep/vm").getSpace().getDimensions();
//     num_nodes = (int)dims_vm[0];

//     auto dims_sv = file.getDataSet("/ep/state_variables").getSpace().getDimensions();
//     num_vars = (int)dims_sv[1];
// }

// void WriterHDF5::read_checkpoint_data(const std::string &filename, double *vm, double *state_vars)
// {
//     HighFive::File file(filename, HighFive::File::ReadOnly);

//     // read_raw fills the caller-provided buffers, which must already be sized
//     // (via read_checkpoint_metadata) to match the datasets on disk.
//     file.getDataSet("/ep/vm").read_raw(vm);
//     file.getDataSet("/ep/state_variables").read_raw(state_vars);
// }

// void WriterHDF5::write_mech_checkpoint(int step, double current_time, int load_increment, double load_factor,
//                                        const double *x_current, const double *fext0, int num_dofs)
// {
//     std::string final_filename = "checkpoint_step_" + std::to_string(step) + ".h5";
//     std::string temp_filename = "checkpoint_mech_tmp.h5";
    
//     hid_t file_id;
//     bool is_new_file = false;

//     file_id = H5Fcreate(temp_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//     is_new_file = true;
    
//     if (is_new_file) {
//         hid_t space_scalar = H5Screate(H5S_SCALAR);
//         hid_t attr_step = H5Acreate(file_id, "step", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//         H5Awrite(attr_step, H5T_NATIVE_INT, &step);
//         H5Aclose(attr_step);

//         hid_t attr_time = H5Acreate(file_id, "time", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//         H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &current_time);
//         H5Aclose(attr_time);
//         H5Sclose(space_scalar);
//     }

//     hid_t group_mech;
//     if (H5Lexists(file_id, "/mechanics", H5P_DEFAULT) > 0) {
//         group_mech = H5Gopen2(file_id, "/mechanics", H5P_DEFAULT);
//     } else {
//         group_mech = H5Gcreate2(file_id, "/mechanics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     }

//     hid_t space_scalar = H5Screate(H5S_SCALAR);
    
//     if (H5Aexists(group_mech, "load_increment")) H5Adelete(group_mech, "load_increment");
//     hid_t attr_load = H5Acreate(group_mech, "load_increment", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_load, H5T_NATIVE_INT, &load_increment);
//     H5Aclose(attr_load);

//     if (H5Aexists(group_mech, "load_factor")) H5Adelete(group_mech, "load_factor");
//     hid_t attr_lf = H5Acreate(group_mech, "load_factor", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_lf, H5T_NATIVE_DOUBLE, &load_factor);
//     H5Aclose(attr_lf);

//     H5Sclose(space_scalar);

//     auto save_array = [&](const char* name, hid_t space, const double* data) {
//         hid_t dset;
//         if (H5Lexists(group_mech, name, H5P_DEFAULT) > 0) {
//             dset = H5Dopen2(group_mech, name, H5P_DEFAULT);
//         } else {
//             dset = H5Dcreate(group_mech, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//         }
//         H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
//         H5Dclose(dset);
//     };

//     hsize_t dims[1] = { (hsize_t)num_dofs };
//     hid_t space_nodal = H5Screate_simple(1, dims, NULL);
    
//     save_array("x_current", space_nodal, x_current);
//     save_array("fext0", space_nodal, fext0);
    
//     H5Sclose(space_nodal);
//     H5Gclose(group_mech);
    
//     H5Fclose(file_id);

//     std::rename(temp_filename.c_str(), final_filename.c_str());
// }

// void WriterHDF5::read_mech_checkpoint_metadata(const std::string &filename, 
//                                                int &step, double &time, int &load_increment, 
//                                                double &load_factor, int &num_dofs)
// {
//     HighFive::File file(filename, HighFive::File::ReadOnly);

//     file.getAttribute("step").read(step);
//     file.getAttribute("time").read(time);

//     auto group_mech = file.getGroup("/mechanics");
//     group_mech.getAttribute("load_increment").read(load_increment);
//     group_mech.getAttribute("load_factor").read(load_factor);

//     auto dims = group_mech.getDataSet("x_current").getSpace().getDimensions();
//     num_dofs = (int)dims[0];
// }

// void WriterHDF5::read_mech_checkpoint_data(const std::string &filename, 
//                                            double *x_current, double *fext0)
// {
//     HighFive::File file(filename, HighFive::File::ReadOnly);
//     auto group_mech = file.getGroup("/mechanics");

//     // read_raw fills the caller-provided buffers, which must already be sized
//     // (via read_mech_checkpoint_metadata) to match the datasets on disk.
//     group_mech.getDataSet("x_current").read_raw(x_current);
//     group_mech.getDataSet("fext0").read_raw(fext0);
// }

// void WriterHDF5::write_coupled_checkpoint(int step, double current_time, 
//                                           const double *vm, const double *state_vars, int num_state_vars,
//                                           int load_increment, double load_factor, 
//                                           const double *x_current, const double *fext0, int num_dofs)
// {
//     char time_buf[64];
//     std::sprintf(time_buf, "%.2f", current_time); 
    
//     std::string final_filename = "checkpoint_t_" + std::string(time_buf) + "ms.h5";
//     std::string temp_filename = "checkpoint_coupled_tmp.h5";

//     hid_t file_id = H5Fcreate(temp_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

//     hid_t space_scalar = H5Screate(H5S_SCALAR);
    
//     hid_t attr_step = H5Acreate(file_id, "step", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_step, H5T_NATIVE_INT, &step);
//     H5Aclose(attr_step);

//     hid_t attr_time = H5Acreate(file_id, "time", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_time, H5T_NATIVE_DOUBLE, &current_time);
//     H5Aclose(attr_time);

    
//     hid_t group_ep = H5Gcreate2(file_id, "/ep", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     hsize_t np = mesh->get_n_points();

//     // vm
//     hsize_t dims_vm[1] = { np };
//     hid_t space_vm = H5Screate_simple(1, dims_vm, NULL);
//     hid_t dataset_vm = H5Dcreate(group_ep, "vm", H5T_NATIVE_DOUBLE, space_vm, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     H5Dwrite(dataset_vm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vm);
//     H5Dclose(dataset_vm);
//     H5Sclose(space_vm);

//     // state_variables
//     hsize_t dims_sv[2] = { np, (hsize_t)num_state_vars };
//     hid_t space_sv = H5Screate_simple(2, dims_sv, NULL);
//     hid_t dataset_sv = H5Dcreate(group_ep, "state_variables", H5T_NATIVE_DOUBLE, space_sv, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     H5Dwrite(dataset_sv, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, state_vars);
//     H5Dclose(dataset_sv);
//     H5Sclose(space_sv);

//     H5Gclose(group_ep);

//     hid_t group_mech = H5Gcreate2(file_id, "/mechanics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

//     hid_t attr_load = H5Acreate(group_mech, "load_increment", H5T_NATIVE_INT, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_load, H5T_NATIVE_INT, &load_increment);
//     H5Aclose(attr_load);

//     hid_t attr_lf = H5Acreate(group_mech, "load_factor", H5T_NATIVE_DOUBLE, space_scalar, H5P_DEFAULT, H5P_DEFAULT);
//     H5Awrite(attr_lf, H5T_NATIVE_DOUBLE, &load_factor);
//     H5Aclose(attr_lf);

//     hsize_t dims_mech[1] = { (hsize_t)num_dofs };
//     hid_t space_mech = H5Screate_simple(1, dims_mech, NULL);

//     hid_t dset_x = H5Dcreate(group_mech, "x_current", H5T_NATIVE_DOUBLE, space_mech, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     H5Dwrite(dset_x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x_current);
//     H5Dclose(dset_x);

//     hid_t dset_fext = H5Dcreate(group_mech, "fext0", H5T_NATIVE_DOUBLE, space_mech, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     H5Dwrite(dset_fext, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fext0);
//     H5Dclose(dset_fext);

//     H5Sclose(space_mech);
//     H5Gclose(group_mech);

//     H5Sclose(space_scalar);
//     H5Fclose(file_id);

//     std::rename(temp_filename.c_str(), final_filename.c_str());
// }