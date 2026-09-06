#include "writer_hdf5.hpp"

WriterHDF5::WriterHDF5(Mesh * m) : Writer(m) //, mesh(m)
{
  // do nothing
} 

WriterHDF5::~WriterHDF5()
{
  close();
}

void WriterHDF5::open(const std::string & file, int nsteps, double step, bool bido, bool is_restart)
{   
    std::size_t pos  = file.find_last_of("/");
    std::string base = file.c_str(); 
     if (pos != std::string::npos) {
        base = file.substr(pos+1);
    }
    if (is_restart) base = base + "_restarted";
    h5name = base + ".h5";

    write_hdf5(file, nsteps, step);
    write_xdmf(file, nsteps, step, bido);
}

void WriterHDF5::close()
{
  // check if need to do something else
}

void WriterHDF5::write_hdf5(const std::string & file, int nsteps, double step)
{
    // prepare to write HDF5 file
    hsize_t np = mesh->get_n_points();
    hsize_t dims[3];
    hid_t file_id, group_id, dataset_id, dataspace_id, props;
    herr_t status;

    std::size_t pos  = file.find_last_of("/");
    std::string base = file.c_str();
    if (pos != std::string::npos)
      base = file.substr(pos+1);

    h5name = base + ".h5";

    double fill_zero = 0.0;

    file_id = H5Fcreate(h5name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // cria o grupo "/geometry" no arquivo
    group_id = H5Gcreate2(file_id, "/geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(group_id);
    
    // cria o grupo "/topology" no arquivo
    group_id = H5Gcreate2(file_id, "/topology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(group_id);

    // cria o grupo "/point_data" no arquivo
    group_id = H5Gcreate2(file_id, "/point_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(group_id);

    // cria o grupo "/cell_data" no arquivo
    group_id = H5Gcreate2(file_id, "/cell_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Gclose(group_id);

    //
    // prepare to write time array to HDF5 file
    //
    double * time = new double[nsteps];
    for(int i=0; i<nsteps; i++) 
        time[i] = i*step;

    // open HDF5 dataset and write array
    file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if(H5Lexists(file_id, "/time", H5P_DEFAULT) == false)
    {
        dims[0] = nsteps;
        dims[1] = 1;
        dims[2] = 1;
        dataspace_id = H5Screate_simple(2, dims, NULL);
       
        dataset_id = H5Dcreate(file_id, "/time", H5T_NATIVE_DOUBLE, dataspace_id,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, time);
        status = H5Dclose(dataset_id);    
        status = H5Sclose(dataspace_id);
    }
    delete time;
      
    //
    // write coordinates of the mesh in HDF5 file
    //
    std::vector<arma::vec3> pts = mesh->get_points();
    double * coords = new double[3*np];
    for(uint i=0; i<np; i++)
    {
        arma::vec3 pt = pts[i];
        coords[i*3 + 0] = pt(0);
        coords[i*3 + 1] = pt(1);
        coords[i*3 + 2] = pt(2);
    }                        
    dims[0] = np;
    dims[1] = 3;
    dataspace_id = H5Screate_simple(2, dims, NULL);   
    dataset_id = H5Dcreate(file_id, "/geometry/coordinates", H5T_NATIVE_DOUBLE, 
                           dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                    H5P_DEFAULT, coords);    
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);
    delete [] coords;     
    
    //
    // write connectivity of the mesh in HDF5 file 
    //
    int ne = mesh->get_n_elements();
    int nn = mesh->get_nen();      
    int * connec = new int[ne*nn];
    for(int i=0; i<ne; i++)
    {   
        std::vector<int> ptnums;
        mesh->get_element_pt_nums(i, ptnums);
        for(int j=0; j<nn; j++)
        {
            connec[i*nn + j] = ptnums[j];
        }
    }

    // TODO: mesh is fixed for one type of element only -> improve this
    dims[0] = ne;    
    dims[1] = nn;
    dataspace_id = H5Screate_simple(2, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/topology/connectivity", H5T_NATIVE_INT, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, connec);
    status = H5Dclose(dataset_id);    
    status = H5Sclose(dataspace_id);
    delete connec;
    
    // 
    // vertex fields 
    //

    // all vertex fields will be initialized to zero
    props = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_fill_value(props, H5T_NATIVE_DOUBLE, &fill_zero);

    // VM array
    //double * vm = new double[nsteps*np];
    //for(hsize_t i=0; i<nsteps*np; i++)
    //    vm[i] = 0.0;
    props = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_fill_value(props, H5T_NATIVE_DOUBLE, &fill_zero);
                
    // cria o dataset coordinates0 (inicial)
    dims[0] = nsteps;
    dims[1] = np;        
    dataspace_id = H5Screate_simple(2, dims, NULL);          
    dataset_id = H5Dcreate(file_id, "/point_data/vm", H5T_NATIVE_DOUBLE, 
                            dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
    //                H5P_DEFAULT, vm);
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);

    // ACTIVE STRESS array (nodal, same layout as vm). The active tension of
    // the cell model lives on the nodes, exactly like the potential.
    dims[0] = nsteps;
    dims[1] = np;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/point_data/active_stress",
                           H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,
                           props, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // FIBRE STRETCH array (nodal, same layout as vm and active_stress).
    // lambda_f = sqrt(I4f) comes from the mechanics and is handed to the
    // cell model; 
    dims[0] = nsteps;
    dims[1] = np;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/point_data/lambda_f",
                           H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,
                           props, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // FIBRE STRETCH RATE array (nodal, same layout as lambda_f).
    // d(lambda_f)/dt in the CELL MODEL's own time unit (1/ms for ToRORd)
    dims[0] = nsteps;
    dims[1] = np;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/point_data/lambda_rate",
                           H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,
                           props, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    // cria o dataset coordinates0 (inicial)
    dims[0] = nsteps;
    dims[1] = ne;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/cell_data/stress", H5T_NATIVE_DOUBLE,
                           dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
    //                H5P_DEFAULT, vm);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/cell_data/strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

    // TENSAO ATIVA EFETIVAMENTE APLICADA (por elemento).
    // Media nodal de Ta no elemento vezes o ta_scale do material. E o que a
    // montagem usa; o campo nodal "active_stress" e o valor BRUTO do modelo
    // celular, antes da escala por material, e os dois nao coincidem.
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/cell_data/Ta_applied",
                           H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,
                           props, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // Mapa do multiplicador ta_scale por elemento. Constante no tempo, mas
    // gravado como campo para poder ser sobreposto aos demais no ParaView.
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/cell_data/ta_scale",
                           H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,
                           props, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/cell_data/long_strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/cell_data/circ_strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/cell_data/rad_strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

  dataspace_id = H5Screate_simple(2, dims, NULL);
  dataset_id = H5Dcreate(file_id, "/cell_data/fibrosis", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
  int * aha = new int[nsteps*ne];
  for(hsize_t i=0; i<nsteps; i++)
    for(hsize_t e=0; e<ne; e++)
      aha[e + ne*i] = mesh->get_element(e).get_index();
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, aha);

  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  //fim teste

  dataspace_id = H5Screate_simple(2, dims, NULL);
  dataset_id = H5Dcreate(file_id, "/cell_data/aha_marker", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

  for(hsize_t i=0; i<nsteps; i++)
    for(hsize_t e=0; e<ne; e++)
      aha[e + ne*i] = mesh->get_element(e).get_aha_num();
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, aha);

  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  //fim teste

  delete [] aha;

    //delete [] vm;
        
    // DISPLACEMENT array
    //double * displ = new double[nsteps*3*np];
    //for(hsize_t i=0; i<nsteps*3*np; i++)
    //    displ[i] = 0.0;

    // cria o dataset coordinates0 (inicial)
    dims[0] = nsteps;
    dims[1] = np;
    dims[2] = 3;
    dataspace_id = H5Screate_simple(3, dims, NULL);       
    dataset_id = H5Dcreate(file_id, "/point_data/displacements", H5T_NATIVE_DOUBLE,
                            dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
    //                H5P_DEFAULT, displ);
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);
    //delete [] displ;
    
    // fecha o arquivo
    status = H5Fclose(file_id);

    if(status != 0) H5Eprint2(status,NULL);
}


void WriterHDF5::write_cell_field_step(int step, const double *data, string fieldname)
{
  hid_t file_id, dataset_id, dataspace_id, memspace_id;
  herr_t status;

  // open an existing file.
  //cout << "Filename: " << h5name.c_str() << endl;
  file_id    = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  string auxstr = string("/cell_data/") + fieldname.c_str();
  dataset_id = H5Dopen(file_id, auxstr.c_str(), H5P_DEFAULT);

  hsize_t k = step;
  hsize_t np = mesh->get_n_elements();
  hsize_t dims[2]   = {1,np};
  hsize_t start[2]  = {k,0};
  hsize_t count[2]  = {1,np};
  hsize_t stride[2] = {1,1};
  hsize_t block[2]  = {1,1};

  // define memory dataspace
  memspace_id = H5Screate_simple(2, dims, NULL);

  // select hyperslab
  dataspace_id = H5Dget_space(dataset_id);
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET,
                               start, stride, count, block);

  // write data
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
                    H5P_DEFAULT, data);

  // close stuff
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  if(status != 0) H5Eprint2(status,NULL);
}

void WriterHDF5::write_eikonal_lat(const std::string & file, const double *lat_data)
{
    std::string base = file;
    std::size_t pos = file.find_last_of("/");
    if (pos != std::string::npos) base = file.substr(pos + 1);
    
    std::string h5_file = base + "_lat.h5";
    std::string xmf_file = base + "_lat.xmf";

    hid_t file_id = H5Fcreate(h5_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Gclose(H5Gcreate2(file_id, "/geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    H5Gclose(H5Gcreate2(file_id, "/topology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    H5Gclose(H5Gcreate2(file_id, "/point_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    hsize_t np = mesh->get_n_points();
    std::vector<arma::vec3> pts = mesh->get_points();
    double *coords = new double[3 * np];
    for (uint i = 0; i < np; i++) {
        coords[i * 3 + 0] = pts[i](0);
        coords[i * 3 + 1] = pts[i](1);
        coords[i * 3 + 2] = pts[i](2);
    }
    hsize_t dims_geom[2] = {np, 3};
    hid_t space_geom = H5Screate_simple(2, dims_geom, NULL);
    hid_t dset_geom = H5Dcreate(file_id, "/geometry/coordinates", H5T_NATIVE_DOUBLE, space_geom, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_geom, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);
    H5Dclose(dset_geom); H5Sclose(space_geom);
    delete[] coords;

    int ne = mesh->get_n_elements();
    int nn = mesh->get_nen();
    int *connec = new int[ne * nn];
    for (int i = 0; i < ne; i++) {
        std::vector<int> ptnums;
        mesh->get_element_pt_nums(i, ptnums);
        for (int j = 0; j < nn; j++) connec[i * nn + j] = ptnums[j];
    }
    hsize_t dims_top[2] = {(hsize_t)ne, (hsize_t)nn};
    hid_t space_top = H5Screate_simple(2, dims_top, NULL);
    hid_t dset_top = H5Dcreate(file_id, "/topology/connectivity", H5T_NATIVE_INT, space_top, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_top, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, connec);
    H5Dclose(dset_top); H5Sclose(space_top);
    delete[] connec;

    hsize_t dims_lat[1] = {np};
    hid_t space_lat = H5Screate_simple(1, dims_lat, NULL);
    hid_t dset_lat = H5Dcreate(file_id, "/point_data/lat", H5T_NATIVE_DOUBLE, space_lat, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_lat, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lat_data);
    H5Dclose(dset_lat); H5Sclose(space_lat);

    H5Fclose(file_id);

    int nd = mesh->get_n_dim();
    std::string toptype;
    if (nn == 2) toptype = "Polyline";
    else if (nn == 3) toptype = "Triangle";
    else if (nn == 4 && nd == 2) toptype = "Quadrilateral";
    else if (nn == 4 && nd == 3) toptype = "Tetrahedron";
    else if (nn == 8) toptype = "Hexahedron";

    std::ofstream xmf(xmf_file.c_str());
    xmf << "<?xml version=\"1.0\" ?>\n"
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
    xmf.close();
}

void WriterHDF5::write_point_field_step(int step, const double *data,
                                        string fieldname)
{
  // Same as write_vm_step, but for any node-centred scalar field. The
  // dataset must already exist in the HDF5 file (see write_hdf5) and the
  // matching attribute must be declared in the XDMF (see write_xdmf).
  hid_t file_id, dataset_id, dataspace_id, memspace_id;
  herr_t status;

  file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  string auxstr = string("/point_data/") + fieldname.c_str();

  // A field that was never created in write_hdf5() cannot be written into.
  // Without this check H5Dopen fails and every call afterwards (get_space,
  // select_hyperslab, write) fails too, burying the real cause under a wall
  // of HDF5-DIAG output once per output step.
  if (H5Lexists(file_id, auxstr.c_str(), H5P_DEFAULT) <= 0)
  {
    cout << " Warning: nodal field '" << fieldname
         << "' has no dataset in the HDF5 file; not written."
         << " Add it to WriterHDF5::write_hdf5() and"
         << " WriterHDF5::write_xdmf()." << endl;
    H5Fclose(file_id);
    return;
  }

  dataset_id = H5Dopen(file_id, auxstr.c_str(), H5P_DEFAULT);

  hsize_t k = step;
  hsize_t np = mesh->get_n_points();
  hsize_t dims[2]   = {1,np};
  hsize_t start[2]  = {k,0};
  hsize_t count[2]  = {1,np};
  hsize_t stride[2] = {1,1};
  hsize_t block[2]  = {1,1};

  // define memory dataspace
  memspace_id = H5Screate_simple(2, dims, NULL);

  // select hyperslab
  dataspace_id = H5Dget_space(dataset_id);
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET,
                               start, stride, count, block);

  // write data
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
                    H5P_DEFAULT, data);

  // close stuff
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  if(status != 0) H5Eprint2(status,NULL);
}


void WriterHDF5::write_vm_step(int step, const double *data)
{    
    hid_t file_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;

    // open an existing file.
    // cout << "VM Filename: " << h5name.c_str() << endl;

    file_id    = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id, "/point_data/vm", H5P_DEFAULT);    
   
    hsize_t k = step;
    hsize_t np = mesh->get_n_points();    
    hsize_t dims[2]   = {1,np};
    hsize_t start[2]  = {k,0};
    hsize_t count[2]  = {1,np};
    hsize_t stride[2] = {1,1};    
    hsize_t block[2]  = {1,1};    
    
    // define memory dataspace
    memspace_id = H5Screate_simple(2, dims, NULL);
  
    // select hyperslab
    dataspace_id = H5Dget_space(dataset_id);   
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, 
            start, stride, count, block);           
    
    // write data
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, 
            H5P_DEFAULT, data);    
     
    // close stuff
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

    if(status != 0) H5Eprint2(status,NULL);
}

void WriterHDF5::write_displ_step(int step, const double *displ)
{
    hid_t file_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;    
       
    // open an existing file.
  //cout <<  "H5Fopen" << endl;
    file_id    = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  //cout <<  "H5Dopen" << endl;
    dataset_id = H5Dopen(file_id, "/point_data/displacements", H5P_DEFAULT);
       
    hsize_t k = step;
    hsize_t np = mesh->get_n_points();    
    hsize_t dims[3]   = {1,np,3};    
    hsize_t start[3]  = {k,0,0};
    hsize_t count[3]  = {1,np,3};
    hsize_t stride[3] = {1,1,1};    
    hsize_t block[3]  = {1,1,1};        
    
    // define memory dataspace
  //cout <<  "H5Screate_simple" << endl;
    memspace_id = H5Screate_simple(3, dims, NULL);
  
    // select hyperslab
  //cout <<  "H5Dget_space" << endl;
    dataspace_id = H5Dget_space(dataset_id);
  //cout <<  "H5Sselect_hyperslab" << endl;
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET,
                                 start, stride, count, block);
  if(status != 0) H5Eprint2(status,NULL);
    // write data
  //cout <<  "H5Dwrite" << endl;
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id,
                      dataspace_id, H5P_DEFAULT, displ);
  if(status != 0) H5Eprint2(status,NULL);
  //cout <<  "H5Sclose" << endl;
  status = H5Sclose(dataspace_id);
  if(status != 0) H5Eprint2(status,NULL);
  //cout <<  "H5Dclose" << endl;
  status = H5Dclose(dataset_id);
  if(status != 0) H5Eprint2(status,NULL);
  //cout <<  "H5Fclose" << endl;
  status = H5Fclose(file_id);

    if(status != 0) H5Eprint2(status,NULL);
}

void WriterHDF5::add_ve()
{       
    hid_t  file_id, dataset_id, dataspace_id;
    herr_t status;

    file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    
    // query vm dataset to find dimensions
    dataset_id = H5Dopen(file_id, "/point_data/vm", H5P_DEFAULT);     
    dataspace_id = H5Dget_space(dataset_id);
    const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);    
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);    

    // then create ve dataset
    hsize_t nsteps = dims[0];
    hsize_t nnodes = dims[1];

    double * ve = new double[nsteps*nnodes];
    for(hsize_t i=0; i<nsteps*nnodes; i++)
        ve[i] = 0.0;
                
    // cria o dataset coordinates0 (inicial)
    dims[0] = nsteps;
    dims[1] = nnodes;        
    dataspace_id = H5Screate_simple(2, dims, NULL);          
    dataset_id = H5Dcreate(file_id, "/point_data/ve", H5T_NATIVE_DOUBLE, 
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                    H5P_DEFAULT, ve);
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

    if(status != 0) H5Eprint2(status,NULL);

    delete [] ve;
}

void WriterHDF5::add_scalar_field(std::string & field_name)
{
    hid_t  file_id, dataset_id, dataspace_id;
    herr_t status;

    file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    std::string scalarfield = "/point_data/" + field_name;

    // query vm dataset to find dimensions
    dataset_id = H5Dopen(file_id, scalarfield.c_str(), H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // then create ve dataset
    hsize_t nsteps = dims[0];
    hsize_t nnodes = dims[1];

    double * ve = new double[nsteps*nnodes];
    for(hsize_t i=0; i<nsteps*nnodes; i++)
        ve[i] = 0.0;

    // cria o dataset coordinates0 (inicial)
    dims[0] = nsteps;
    dims[1] = nnodes;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/point_data/ve", H5T_NATIVE_DOUBLE,
                           dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, ve);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

    if(status != 0) H5Eprint2(status,NULL);

    delete [] ve;
}

void WriterHDF5::add_fibers()
{
    /*
    hid_t  file_id, dataset_id, dataspace_id;
    herr_t status;

    file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    // query vm dataset to find dimensions
    dataset_id = H5Dopen(file_id, "/point_data/vm", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    // then create ve dataset
    hsize_t nsteps = dims[0];
    hsize_t nnodes = dims[1];
    double * ve = new double[nsteps*nnodes];
    for(hsize_t i=0; i<nsteps*nnodes; i++)
        ve[i] = 0.0;

    // cria o dataset coordinates0 (inicial)
    dims[0] = nsteps;
    dims[1] = nnodes;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/point_data/ve", H5T_NATIVE_DOUBLE,
                           dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, ve);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

    delete [] ve;
     */
}

void WriterHDF5::write_ve_step(int step, const double *data)
{    
    hid_t     file_id, dataset_id, dataspace_id, memspace_id;
    herr_t    status;    

    // open an existing file.
    file_id    = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id, "/point_data/ve", H5P_DEFAULT);    
   
    hsize_t k = step;
    hsize_t np = mesh->get_n_points();    
    hsize_t dims[2]   = {1,np};
    hsize_t start[2]  = {k,0};
    hsize_t count[2]  = {1,np};
    hsize_t stride[2] = {1,1};    
    hsize_t block[2]  = {1,1};    
    
    // define memory dataspace
    memspace_id = H5Screate_simple(2, dims, NULL);
  
    // select hyperslab
    dataspace_id = H5Dget_space(dataset_id);   
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, 
            start, stride, count, block);           
    
    // write data
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, 
            H5P_DEFAULT, data);    
     
    // close stuff
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);

    if(status != 0) H5Eprint2(status,NULL);
}

void WriterHDF5::write_xdmf(const std::string & file, int nsteps, 
                            double step, bool bido)
{
    int np = mesh->get_n_points();
    int ne = mesh->get_n_elements();
    int nn = mesh->get_nen(); 
    int nd = mesh->get_n_dim();

    std::string toptype;
    if (nn == 2) toptype = "Polyline";
    else if (nn == 3) toptype = "Triangle";
    else if (nn == 4 && nd == 2) toptype = "Quadrilateral";
    else if (nn == 4 && nd == 3) toptype = "Tetrahedron";
    else if (nn == 8) toptype = "Hexahedron";

    std::size_t pos  = file.find_last_of("/");
    std::string base = file.c_str();
    if (pos != std::string::npos)
      base = file.substr(pos+1);
    std::string xmf_file(base + ".xmf");
    
    if(bido)
    {
        // do nothing for now...
    }
        
    //
    // open the file and write the XML description of the mesh
    //
    xmf.open(xmf_file.c_str());
    xmf << "<?xml version=\"1.0\" ?>\n"
        << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
        << "<Xdmf Version=\"2.0\">\n"
        << "  <Domain>\n\n";
            
    xmf << "    <Topology TopologyType=\"" << toptype << "\"\n"
        << "        NumberOfElements=\"" << ne << "\">\n"
        << "    <DataItem Format=\"HDF\" \n"
        << "        DataType=\"Int\"  \n"
        << "        Dimensions=\"" << ne << " " << nn << "\">" 
        << h5name << ":/topology/connectivity\n"
        << "    </DataItem>\n"
        << "    </Topology>\n\n";        
    
    xmf << "    <Geometry GeometryType=\"XYZ\">\n"
        << "        <DataItem Dimensions=\"" << np << " " << 3 << "\"\n"
        << "            NumberType=\"Double\" \n"
        << "            Precision=\"8\" \n"
        << "            Format=\"HDF\">" << h5name << ":/geometry/coordinates\n"
        << "        </DataItem>\n"
        << "    </Geometry>\n\n";        
    
    //
    // write time information
    //    
    xmf << "    <Grid Name=\"TimeSeries\" \n"
        << "        GridType=\"Collection\"\n"
        << "        CollectionType=\"Temporal\">\n"           
        << "        <Time TimeType=\"List\">\n" // for NON-UNIFORM time steps
        << "            <DataItem  Format=\"HDF\" \n" 
        << "                NumberType=\"Double\" \n"
        << "                Dimensions=\"" << nsteps << "\">" << h5name << ":/time\n"
        << "            </DataItem>\n"
        << "        </Time>\n";
    
    for(int i=0; i<nsteps; i++)
    {   
        xmf << "        <Grid Name=\"T" << i << "\" GridType=\"Uniform\">\n"
            << "            <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n"
            << "            <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n"; 
        
        //
        // potential
        //     
        xmf << "            <Attribute Name=\"vm\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Node\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << np << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << np <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Points\" \n"
            << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/point_data/vm\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";

        //
        // active tension (nodal)
        //
        xmf << "            <Attribute Name=\"active_stress\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Node\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << np << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << np <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Points\" \n"
            << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/point_data/active_stress\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";

        //
        // fibre stretch lambda_f (nodal)
        //
        xmf << "            <Attribute Name=\"lambda_f\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Node\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << np << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << np <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Points\" \n"
            << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/point_data/lambda_f\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";

        //
        // fibre stretch rate d(lambda_f)/dt (nodal)
        //
        xmf << "            <Attribute Name=\"lambda_rate\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Node\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << np << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << np <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Points\" \n"
            << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/point_data/lambda_rate\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";

        //teste cell field
        //Ta_applied -- tensao ativa efetivamente aplicada pela mecanica
        xmf << "            <Attribute Name=\"Ta_applied\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Cell\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << ne << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << ne <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Cells\" \n"
            << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/cell_data/Ta_applied\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        //ta_scale -- mapa do multiplicador por elemento
        xmf << "            <Attribute Name=\"ta_scale\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Cell\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << ne << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << ne <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Cells\" \n"
            << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/cell_data/ta_scale\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        //stress
        xmf << "            <Attribute Name=\"stress\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Cell\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << ne << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << ne <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Cells\" \n"
            << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/cell_data/stress\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        //strain
        xmf << "            <Attribute Name=\"strain\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Cell\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << ne << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << ne <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Cells\" \n"
            << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/cell_data/strain\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        //long_strain
        xmf << "            <Attribute Name=\"long_strain\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Cell\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << ne << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << ne <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Cells\" \n"
            << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/cell_data/long_strain\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        //circ_strain
        xmf << "            <Attribute Name=\"circ_strain\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Cell\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << ne << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << ne <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Cells\" \n"
            << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/cell_data/circ_strain\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        //rad_strain
        xmf << "            <Attribute Name=\"rad_strain\" \n"
            << "                AttributeType=\"Scalar\" \n"
            << "                Center=\"Cell\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << ne << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
            << "                    " << i << " 0 \n"
            << "                    1 1 \n"
            << "                    1 " << ne <<"\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Cells\" \n"
            << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/cell_data/rad_strain\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
  	//fibrosis
      xmf << "            <Attribute Name=\"fibrosis\" \n"
          << "                AttributeType=\"Scalar\" \n"
          << "                Center=\"Cell\">\n"
          << "            <DataItem ItemType=\"HyperSlab\" \n"
          << "                Dimensions=\"1 " << ne << "\" \n"
          << "                Type=\"HyperSlab\">\n"
          << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
          << "                    " << i << " 0 \n"
          << "                    1 1 \n"
          << "                    1 " << ne <<"\n"
          << "                </DataItem>\n"
          << "                <DataItem Name=\"Cells\" \n"
          << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
          << "                    Format=\"HDF\">" << h5name << ":/cell_data/fibrosis\n"
          << "                </DataItem>\n"
          << "            </DataItem>\n"
          << "            </Attribute>\n";
    //aha_segment
      xmf << "            <Attribute Name=\"aha_marker\" \n"
          << "                AttributeType=\"Scalar\" \n"
          << "                Center=\"Cell\">\n"
          << "            <DataItem ItemType=\"HyperSlab\" \n"
          << "                Dimensions=\"1 " << ne << "\" \n"
          << "                Type=\"HyperSlab\">\n"
          << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
          << "                    " << i << " 0 \n"
          << "                    1 1 \n"
          << "                    1 " << ne <<"\n"
          << "                </DataItem>\n"
          << "                <DataItem Name=\"Cells\" \n"
          << "                    Dimensions=\"" << nsteps << " " << ne << "\" \n"
          << "                    Format=\"HDF\">" << h5name << ":/cell_data/aha_marker\n"
          << "                </DataItem>\n"
          << "            </DataItem>\n"
          << "            </Attribute>\n";
          //fim teste
        
        if(bido)
        {
            //
            // extracellular potential
            //     
            xmf << "            <Attribute Name=\"ve\" \n"
                << "                AttributeType=\"Scalar\" \n"
                << "                Center=\"Node\">\n"
                << "            <DataItem ItemType=\"HyperSlab\" \n"
                << "                Dimensions=\"1 " << np << "\" \n"
                << "                Type=\"HyperSlab\">\n"
                << "                <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
                << "                    " << i << " 0 \n"
                << "                    1 1 \n"
                << "                    1 " << np <<"\n"
                << "                </DataItem>\n"
                << "                <DataItem Name=\"Points\" \n"
                << "                    Dimensions=\"" << nsteps << " " << np << "\" \n"
                << "                    Format=\"HDF\">" << h5name << ":/point_data/ve\n"
                << "                </DataItem>\n"
                << "            </DataItem>\n"
                << "            </Attribute>\n";          
        }
        
        //
        // displacements
        //
        xmf << "            <Attribute Name=\"displacement\" \n"
            << "                AttributeType=\"Vector\" \n"
            << "                Center=\"Node\">\n"
            << "            <DataItem ItemType=\"HyperSlab\" \n"
            << "                Dimensions=\"1 " << np  << " " << 3 << "\" \n"
            << "                Type=\"HyperSlab\">\n"
            << "                <DataItem Dimensions=\"3 3\" Format=\"XML\">\n"
            << "                    " << i << " 0 0 \n"
            << "                    1 1 1 \n"
            << "                    1 " << np << " 3" << "\n"
            << "                </DataItem>\n"
            << "                <DataItem Name=\"Points\" \n"
            << "                    Dimensions=\"" << nsteps << " " << np << " 3" << "\" \n"
            << "                    Format=\"HDF\">" << h5name << ":/point_data/displacements\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        xmf << "        </Grid>\n";
    }
    
    // Used to write only the mesh
    //fprintf(xmf, "    <Grid Name=\"Mesh\" GridType=\"Uniform\">\n");
    //fprintf(xmf, "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n");
    //fprintf(xmf, "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");
    //fprintf(xmf, "    </Grid>\n\n");     
    
    //
    // close file
    //
    xmf << "    </Grid>\n"
        << "  </Domain>\n"
        << "</Xdmf>\n";
    
    xmf.close();
}

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

    // saving vm
    hsize_t dims_vm[1] = { np };
    hid_t dataspace_vm = H5Screate_simple(1, dims_vm, NULL);
    hid_t dataset_vm = H5Dcreate(group_ep, "vm", H5T_NATIVE_DOUBLE, dataspace_vm, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dataset_vm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vm);
    
    H5Dclose(dataset_vm);
    H5Sclose(dataspace_vm);

    // saving state_variables
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
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        throw std::runtime_error("Error: Not possible to open checkpoint file: " + filename);
    }

    hid_t attr_step = H5Aopen(file_id, "step", H5P_DEFAULT);
    H5Aread(attr_step, H5T_NATIVE_INT, &step);
    H5Aclose(attr_step);

    hid_t attr_time = H5Aopen(file_id, "time", H5P_DEFAULT);
    H5Aread(attr_time, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attr_time);

    hid_t dataset_vm = H5Dopen2(file_id, "/ep/vm", H5P_DEFAULT);
    hid_t space_vm = H5Dget_space(dataset_vm);
    hsize_t dims_vm[1];
    H5Sget_simple_extent_dims(space_vm, dims_vm, NULL);
    num_nodes = (int)dims_vm[0];
    
    H5Sclose(space_vm);
    H5Dclose(dataset_vm);

    hid_t dataset_sv = H5Dopen2(file_id, "/ep/state_variables", H5P_DEFAULT);
    hid_t space_sv = H5Dget_space(dataset_sv);
    hsize_t dims_sv[2];
    H5Sget_simple_extent_dims(space_sv, dims_sv, NULL);
    num_vars = (int)dims_sv[1];
    
    H5Sclose(space_sv);
    H5Dclose(dataset_sv);
    H5Fclose(file_id);
}

void WriterHDF5::read_checkpoint_data(const std::string &filename, double *vm, double *state_vars)
{
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t dataset_vm = H5Dopen2(file_id, "/ep/vm", H5P_DEFAULT);
    H5Dread(dataset_vm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vm);
    H5Dclose(dataset_vm);

    hid_t dataset_sv = H5Dopen2(file_id, "/ep/state_variables", H5P_DEFAULT);
    H5Dread(dataset_sv, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, state_vars);
    H5Dclose(dataset_sv);

    H5Fclose(file_id);
}

void WriterHDF5::write_mech_checkpoint(int step, double current_time, int load_increment, double load_factor,
                                       const double *x_current, const double *fext0, int num_dofs)
{
    std::string final_filename = "checkpoint_step_" + std::to_string(step) + ".h5";
    std::string temp_filename = "checkpoint_mech_tmp.h5";
    
    hid_t file_id;
    bool is_new_file = false;

    file_id = H5Fcreate(temp_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    is_new_file = true;
    
    if (is_new_file) {
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
    if (H5Lexists(file_id, "/mechanics", H5P_DEFAULT) > 0) {
        group_mech = H5Gopen2(file_id, "/mechanics", H5P_DEFAULT);
    } else {
        group_mech = H5Gcreate2(file_id, "/mechanics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

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
        if (H5Lexists(group_mech, name, H5P_DEFAULT) > 0) {
            dset = H5Dopen2(group_mech, name, H5P_DEFAULT);
        } else {
            dset = H5Dcreate(group_mech, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }
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
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    hid_t attr_step = H5Aopen(file_id, "step", H5P_DEFAULT);
    H5Aread(attr_step, H5T_NATIVE_INT, &step);
    H5Aclose(attr_step);

    hid_t attr_time = H5Aopen(file_id, "time", H5P_DEFAULT);
    H5Aread(attr_time, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attr_time);

    hid_t group_mech = H5Gopen2(file_id, "/mechanics", H5P_DEFAULT);
    
    hid_t attr_load = H5Aopen(group_mech, "load_increment", H5P_DEFAULT);
    H5Aread(attr_load, H5T_NATIVE_INT, &load_increment);
    H5Aclose(attr_load);

    hid_t attr_lf = H5Aopen(group_mech, "load_factor", H5P_DEFAULT);
    H5Aread(attr_lf, H5T_NATIVE_DOUBLE, &load_factor);
    H5Aclose(attr_lf);

    hid_t dset_x = H5Dopen2(group_mech, "x_current", H5P_DEFAULT);
    hid_t space_x = H5Dget_space(dset_x);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(space_x, dims, NULL);
    num_dofs = (int)dims[0];
    
    H5Sclose(space_x);
    H5Dclose(dset_x);
    H5Gclose(group_mech);
    H5Fclose(file_id);
}

void WriterHDF5::read_mech_checkpoint_data(const std::string &filename, 
                                           double *x_current, double *fext0)
{
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group_mech = H5Gopen2(file_id, "/mechanics", H5P_DEFAULT);

    auto read_dataset = [&](const char* name, double* data) {
        hid_t dset = H5Dopen2(group_mech, name, H5P_DEFAULT);
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(dset);
    };

    read_dataset("x_current", x_current);
    read_dataset("fext0", fext0);

    H5Gclose(group_mech);
    H5Fclose(file_id);
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

    // vm
    hsize_t dims_vm[1] = { np };
    hid_t space_vm = H5Screate_simple(1, dims_vm, NULL);
    hid_t dataset_vm = H5Dcreate(group_ep, "vm", H5T_NATIVE_DOUBLE, space_vm, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_vm, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vm);
    H5Dclose(dataset_vm);
    H5Sclose(space_vm);

    // state_variables
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