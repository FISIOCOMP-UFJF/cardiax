#include "writer_hdf5.hpp"

WriterHDF5::WriterHDF5(Mesh * m) : Writer(m) //, mesh(m)
{
  // do nothing
} 

WriterHDF5::~WriterHDF5()
{
  close();
}

void WriterHDF5::open(const std::string & file, int nsteps, double step, bool bido)
{

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

    // cria o grupo "/vertex_field" no arquivo
    group_id = H5Gcreate2(file_id, "/vertex_field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    dataset_id = H5Dcreate(file_id, "/vertex_field/vm", H5T_NATIVE_DOUBLE, 
                            dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
    //                H5P_DEFAULT, vm);
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);

    //teste cell field
    // cria o dataset coordinates0 (inicial)
    dims[0] = nsteps;
    dims[1] = ne;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/vertex_field/stress", H5T_NATIVE_DOUBLE,
                           dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
    //                H5P_DEFAULT, vm);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/vertex_field/strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/vertex_field/long_strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/vertex_field/circ_strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/vertex_field/rad_strain", H5T_NATIVE_DOUBLE,
                         dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    //fim teste

  dataspace_id = H5Screate_simple(2, dims, NULL);
  dataset_id = H5Dcreate(file_id, "/vertex_field/aha_marker", H5T_NATIVE_DOUBLE,
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
    dataset_id = H5Dcreate(file_id, "/vertex_field/displacements", H5T_NATIVE_DOUBLE,
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
  string auxstr = string("/vertex_field/") + fieldname.c_str();
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

void WriterHDF5::write_vm_step(int step, const double *data)
{    
    hid_t file_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;

    // open an existing file.
  cout << "VM Filename: " << h5name.c_str() << endl;

    file_id    = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id, "/vertex_field/vm", H5P_DEFAULT);    
   
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
    dataset_id = H5Dopen(file_id, "/vertex_field/displacements", H5P_DEFAULT);
       
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
    dataset_id = H5Dopen(file_id, "/vertex_field/vm", H5P_DEFAULT);     
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
    dataset_id = H5Dcreate(file_id, "/vertex_field/ve", H5T_NATIVE_DOUBLE, 
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

    std::string scalarfield = "/vertex_field/" + field_name;

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
    dataset_id = H5Dcreate(file_id, "/vertex_field/ve", H5T_NATIVE_DOUBLE,
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
    dataset_id = H5Dopen(file_id, "/vertex_field/vm", H5P_DEFAULT);
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
    dataset_id = H5Dcreate(file_id, "/vertex_field/ve", H5T_NATIVE_DOUBLE,
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
    dataset_id = H5Dopen(file_id, "/vertex_field/ve", H5P_DEFAULT);    
   
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/vm\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";

        //teste cell field
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/stress\n"
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/strain\n"
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/long_strain\n"
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/circ_strain\n"
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/rad_strain\n"
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
          << "                    Format=\"HDF\">" << h5name << ":/vertex_field/aha_marker\n"
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
                << "                    Format=\"HDF\">" << h5name << ":/vertex_field/ve\n"
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/displacements\n"
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

// BACKUP ----------------------------------------------------------------------

/*
void XdmfIO::write_xdmf(const std::string & file, int nsteps, bool mov)
{
  
    xmf = fopen(xmf_file.c_str(), "w");    
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, "  <Domain>\n\n");
            
    fprintf(xmf, "    <Topology TopologyType=\"%s\" 
                          NumberOfElements=\"%d\" >\n", toptype.c_str(), ne);
    fprintf(xmf, "      <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"%d %d\">%s.h5:/topology/connectivity</DataItem>\n", ne, nn, file.c_str());
    fprintf(xmf, "    </Topology>\n\n");
    
    fprintf(xmf, "    <Geometry GeometryType=\"XYZ\">\n");
    fprintf(xmf, "      <DataItem Dimensions=\"%d %d\" NumberType=\"Double\" Precision=\"8\" Format=\"HDF\">%s.h5:/geometry/coordinates</DataItem>\n", np, 3, file.c_str());
    fprintf(xmf, "    </Geometry>\n\n");    
    
    // Used to write only the mesh
    //fprintf(xmf, "    <Grid Name=\"Mesh\" GridType=\"Uniform\">\n");
    //fprintf(xmf, "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n");
    //fprintf(xmf, "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");
    //fprintf(xmf, "    </Grid>\n\n");    
    
    // moved to add_time()
    fprintf(xmf, "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
    fprintf(xmf, "      <Time TimeType=\"List\">\n");
    fprintf(xmf, "        <DataItem  Format=\"HDF\" NumberType=\"Double\" Dimensions=\"%d\">%s.h5:/time</DataItem>\n", NSTEPS, file.c_str());
    fprintf(xmf, "      </Time>\n");      
    
    for(int i=0; i<NSTEPS; i++)
    {   
        fprintf(xmf, "      <Grid Name=\"T%d\" GridType=\"Uniform\">\n", i);
        fprintf(xmf, "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n");
        fprintf(xmf, "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");        
        fprintf(xmf, "        <Attribute Name=\"solution\" Center=\"Node\">\n");
        fprintf(xmf, "          <DataItem Format=\"HDF\" Dimensions=\"%d\">%s.h5:/data/step%d</DataItem>\n", np, file.c_str(), i);
        fprintf(xmf, "        </Attribute>\n");                       
        fprintf(xmf, "      </Grid>\n");
    } 
    
    // moved to close())
    fprintf(xmf, "    </Grid>\n");
    fprintf(xmf, "  </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);    
}

 
void XdmfIO::write_hdf5_displ_step(int step, double * displ)
{
    
    hid_t file_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;
    hsize_t np = mesh->get_n_points();      
    hsize_t dims[2] = {np,3};   
    
    // Open an existing file
    char nome[256];
    
    file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    sprintf(nome,"/displacements/step%d",step);    
    dims[0] = mesh->get_n_points(); 
    dims[1] = 3;
    
    dataspace_id = H5Screate_simple(2, dims, NULL);    
    dataset_id   = H5Dcreate2(file_id, nome, H5T_NATIVE_DOUBLE, dataspace_id, 
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
            H5P_DEFAULT, displ);    
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
}

 
void XdmfIO::write_hdf5_data_step(double * data, int step)
{
    hid_t     file_id, dataset_id, dataspace_id;
    hsize_t   dims[2];
    herr_t    status;
    char      nome[256];

    // Open an existing file.
    file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    sprintf(nome,"/vertex_field/step%d",step);    
    dims[0] = mesh->get_n_points(); 
    dims[1] = 1;
    
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id   = H5Dcreate2(file_id, nome, H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
}

 
void XdmfIO::write_hdf5_displ_step(int step, double * displ)
{
    hid_t file_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;    
       
    // open an existing file.
    file_id    = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id, "/vertex_field/displacements", H5P_DEFAULT);      
       
    hsize_t k = step;
    hsize_t np = mesh->get_n_points();    
    hsize_t dims[3]   = {1,np,3};    
    hsize_t start[3]  = {k,0,0};
    hsize_t count[3]  = {1,np,3};
    hsize_t stride[3] = {1,1,1};    
    hsize_t block[3]  = {1,1,1};    
 
    // Open an existing file.
    char nome[256];
    file_id = H5Fopen(h5name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    sprintf(nome,"/data/displ%d",step);    
    dims[0] = mesh->get_n_points(); 
    dims[1] = 3;
    
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id   = H5Dcreate2(file_id, nome, H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, displ);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    status = H5Fclose(file_id);
}
  
 
void add_time()
{
    fprintf(xmf, "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
    fprintf(xmf, "      <Time TimeType=\"List\">\n");
    fprintf(xmf, "        <DataItem  Format=\"HDF\" NumberType=\"Double\" Dimensions=\"%d\">%s:/time</DataItem>\n", nsteps, h5name.c_str());
    fprintf(xmf, "      </Time>\n");          
    for(int i=0; i<nsteps; i++)
    {   
        fprintf(xmf, "      <Grid Name=\"T%d\" GridType=\"Uniform\">\n", i);
        fprintf(xmf, "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n");
        fprintf(xmf, "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n");                
        
        //
        // potential
        //
        //fprintf(xmf, "        <Attribute Name=\"solution\" Center=\"Node\">\n");
        //fprintf(xmf, "          <DataItem Format=\"HDF\" Dimensions=\"%d\">%s:/data/step%d</DataItem>\n", np, h5name.c_str(), i);
        //fprintf(xmf, "        </Attribute>\n");        
        fprintf(xmf, "<Attribute Name=\"vm\" AttributeType=\"Scalar\" Center=\"Node\">\n");    
        fprintf(xmf, "<DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d\" Type=\"HyperSlab\">\n", 1, np);
        fprintf(xmf, "<DataItem Dimensions=\"3 2\" Format=\"XML\">\n");
        fprintf(xmf, " %d 0 \n", i);
        fprintf(xmf, " 1 1 \n");
        fprintf(xmf, " 1 %d\n", np);
        fprintf(xmf, "</DataItem>\n");
        fprintf(xmf, "<DataItem Name=\"Points\" Dimensions=\"%d %d\" 
                       Format=\"HDF\">%s:/vertex_field/potential</DataItem>\n", nsteps, np, h5name.c_str());
        fprintf(xmf, "</DataItem>\n");        
        fprintf(xmf, "</Attribute>\n");
        
        //
        // displacements
        //        
        fprintf(xmf, "<Attribute Name=\"displacement\" AttributeType=\"Vector\" Center=\"Node\">\n");    
        fprintf(xmf, "<DataItem ItemType=\"HyperSlab\" Dimensions=\"%d %d %d\" Type=\"HyperSlab\">\n", 1, np, 3);
        fprintf(xmf, "<DataItem Dimensions=\"3 3\" Format=\"XML\">\n");
        fprintf(xmf, " %d 0 0 \n", i);
        fprintf(xmf, " 1 1 1 \n");
        fprintf(xmf, " 1 %d 3 \n", np);
        fprintf(xmf, "</DataItem>\n");
        fprintf(xmf, "<DataItem Name=\"Points\" Dimensions=\"%d %d %d\" 
                        Format=\"HDF\">%s:/vertex_field/displacements</DataItem>\n", nsteps, np, 3, h5name.c_str());
        fprintf(xmf, "</DataItem>\n");
        fprintf(xmf, "</Attribute>\n");
        
        fprintf(xmf, "      </Grid>\n");
    } 
 }    
  
*/

/*
void WriterHDF5::write_time(int nsteps, double step)
{    
    // prepare XDMF file structure
    int np = mesh->get_n_points(); 
    
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/vm\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
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
            << "                    Format=\"HDF\">" << h5name << ":/vertex_field/displacements\n"
            << "                </DataItem>\n"
            << "            </DataItem>\n"
            << "            </Attribute>\n";
        xmf << "        </Grid>\n";
    }        
}
*/
