#include "writer_vtk.hpp"

WriterVTK::WriterVTK(Mesh * m): Writer(m)
{
  // do nothing
}

WriterVTK::~WriterVTK()
{
  // do nothing
}

void WriterVTK::open(const std::string & file, int nsteps, double step, bool bido)
{
    //write_hdf5(file, nsteps, step);
    write_xml(file, nsteps, step, bido);
}

void WriterVTK::close()
{

}

void WriterVTK::write_vm_step(int step, const double *data)
{

}

void WriterVTK::write_displ_step(int step, const double *data)
{

}

void WriterVTK::write_xml(const std::string & file, int nsteps,
                          double step, bool bido)
{
  int np = mesh->get_n_points();
  int ne = mesh->get_n_elements();
  int nn = mesh->get_nen();
  int nd = mesh->get_n_dim();
  int type;

  switch(nn)
  {
    case 2: type = 1; break;
    case 3: type = 5; break;
    case 4:
      if(nd == 2) type = 9;
      else if(nd==3) type = 10;
      break;
    case 8: type = 11; break;
    default:
      error("element not found");
  }

  std::string filename;
  filename = file + ".vtu";

  int offset=0;

  if(bido) { }

  xmlfile.open(filename.c_str(), ios::out | ios::binary);
  xmlfile << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          << "  <UnstructuredGrid>\n"
          << "    <Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"" << ne << "\" >\n";
  /*
  xmlfile << "      <PointData Scalars=\"scalars\">\n";
  xmlfile << "      </PointData>\n";

  xmlfile << "      <CellData Scalars=\"scalars\">\n";
  xmlfile << "      </CellData>\n";
  */

  xmlfile << "      <Points>\n"
          << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" \n"
          << "                   format=\"appended\" offset=\"" << offset << "\" />\n";
          //<< "        </DataArray>\n";
  xmlfile << "      </Points>\n";

  offset += 3*np*sizeof(double);

  xmlfile << "      <Cells>\n"
          << "        <DataArray type=\"Int32\" Name=\"connectivity\" \n"
          << "                   format=\"appended\" offset=\"" << offset << "\" />\n";
            // ne linhas por nn ints
          //<< "        </DataArray>\n"

  offset += nn*ne*sizeof(int);

  xmlfile << "        <DataArray type=\"Int32\" Name=\"offsets\" \n"
          << "                   format=\"appended\" offset=\"" << offset << "\" />\n";
          //8 16 24 32 40 48 56 64
          //<< "        </DataArray>\n"

  offset += ne*sizeof(int);

  xmlfile << "        <DataArray type=\"UInt8\" Name=\"types\" \n"
          << "                   format=\"appended\" offset=\"" << offset << "\" />\n"
          // 11 11 11 11 11 11 11 11
          //<< "        </DataArray>\n"
          << "      </Cells>\n";

  offset += ne*sizeof(int);

  xmlfile << "    </Piece>\n";
  xmlfile << "  </UnstructuredGrid>\n";

  xmlfile << "  <AppendedData encoding=\"raw\">\n";
  xmlfile.close();

  write_bin(filename);

  xmlfile.open(filename.c_str(), ios::app | ios::binary);
  xmlfile << "\n";
  xmlfile << "  </AppendedData>\n";
  xmlfile << "</VTKFile>";
  xmlfile.close();
}

void WriterVTK::write_bin(const std::string & filename)
{
  // write coordinates of the mesh
  int np = mesh->get_n_points();
  std::vector<arma::vec3> pts = mesh->get_points();
  double * coords = new double[3*np];
  for(uint i=0; i<np; i++)
  {
    arma::vec3 pt = pts[i];
    coords[i*3 + 0] = pt(0);
    coords[i*3 + 1] = pt(1);
    coords[i*3 + 2] = pt(2);
  }

  // write connectivity of the mesh in HDF5 file
  int ne = mesh->get_n_elements();
  int nn = mesh->get_nen();
  int nd = mesh->get_n_dim();

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

  int type;
  switch(nn)
  {
    case 2: type = 1; break;
    case 3: type = 5; break;
    case 4:
      if(nd == 2) type = 9;
      else if(nd==3) type = 10;
          break;
    case 8: type = 11; break;
    default:
      error("element not found");
  }

  char start = '_';
  int * vtype = new int[ne];
  int * voffset = new int[ne];
  for(int i=0; i<ne; i++) vtype[i] = type;
  cout << "TYPE=" << type << endl;
  for(int i=0; i<ne; i++)
  {
    voffset[i] = (i + 1) * nn;
    cout << "voff_" << i << " " << voffset[i] << endl;
  }

  xmlfile.open(filename.c_str(), ios::app | ios::binary);
  xmlfile.write((char *) &start, sizeof(char));
  xmlfile.write((char *) &coords, 3*np*sizeof(double));
  xmlfile.write((char *) &connec, nn*ne*sizeof(int));
  xmlfile.write((char *) &voffset, ne*sizeof(int));
  xmlfile.write((char *) &vtype, ne*sizeof(int));
  xmlfile.close();

  delete [] vtype;
  delete [] coords;
  delete [] connec;
}

/*
void WriterVTK::set_points(const std::vector<arma::vec3> & pts)
{
  arma::vec3 pt;
  vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();

  for(uint i=0; i<pts.size(); i++)
  {
    pt = pts[i];
    points->InsertNextPoint(pt(0),pt(1),pt(2));
  }
  unstructuredGrid->SetPoints(points);
}

void WriterVTK::set_cells()
{
  int n_dim  = mesh.get_ndim();
  int n_node = mesh.get_nen();
  int n_elem = mesh.get_n_elements();

  std::vector<int> ptnums;

  vtkSmartPointer<vtkCellArray> cellArray =
      vtkSmartPointer<vtkCellArray>::New();

  if (n_dim == 2 && n_node == 3)
  {
    vtkSmartPointer<vtkTriangle> tri =
        vtkSmartPointer<vtkTriangle>::New();

    for(int i=0; i<n_elem; i++)
    {
      mesh.get_element_pt_nums(i, ptnums);
      for(uint k=0; k<ptnums.size(); k++)
        tri->GetPointIds()->SetId(k,ptnums[k]);

      cellArray->InsertNextCell(tri);
    }
    unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);
  }
  else if (n_dim == 2 && n_node == 4)
  {
    vtkSmartPointer<vtkQuad> quad =
        vtkSmartPointer<vtkQuad>::New();

    for(int i=0; i<n_elem; i++)
    {
      mesh.get_element_pt_nums(i, ptnums);
      for(uint k=0; k<ptnums.size(); k++)
        quad->GetPointIds()->SetId(k,ptnums[k]);
      cellArray->InsertNextCell(quad);
    }
    unstructuredGrid->SetCells(VTK_QUAD, cellArray);
  }
  else if (n_dim == 3 && n_node == 8)
  {
    vtkSmartPointer<vtkHexahedron> hex =
        vtkSmartPointer<vtkHexahedron>::New();

    for(int i=0; i<n_elem; i++)
    {
      mesh.get_element_pt_nums(i, ptnums);
      for(uint k=0; k<ptnums.size(); k++)
        hex->GetPointIds()->SetId(k,ptnums[k]);
      cellArray->InsertNextCell(hex);
    }
    unstructuredGrid->SetCells(VTK_HEXAHEDRON, cellArray);
  }
  else if (n_dim == 3 && n_node == 4)
  {
    vtkSmartPointer<vtkTetra> tet =
        vtkSmartPointer<vtkTetra>::New();

    for(int i=0; i<n_elem; i++)
    {
      mesh.get_element_pt_nums(i, ptnums);
      for(uint k=0; k<ptnums.size(); k++)
        tet->GetPointIds()->SetId(k,ptnums[k]);
      cellArray->InsertNextCell(tet);
    }
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);
  } 
  else
  {
    throw std::runtime_error("Invalid element type in VtkIO.");
  }
}

void WriterVTK::set_cells(const std::vector< std::vector<uint> > & connec)
{
  int n_dim = mesh.get_ndim();
  int n_node = mesh.get_nen();

  std::vector<int> ptnums;

  vtkSmartPointer<vtkCellArray> cellArray =
    vtkSmartPointer<vtkCellArray>::New();

  if (n_dim == 2 && n_node == 3)
  {
    vtkSmartPointer<vtkTriangle> tri = 
      vtkSmartPointer<vtkTriangle>::New();

    for(uint i=0; i<connec.size(); i++)
    {
      for(uint k=0; k<connec[i].size(); k++)
        tri->GetPointIds()->SetId(k,connec[i][k]);

      cellArray->InsertNextCell(tri);
    }
    unstructuredGrid->SetCells(VTK_TRIANGLE, cellArray);

  }
  else if (n_dim == 2 && n_node == 4)
  {
    vtkSmartPointer<vtkQuad> quad = 
      vtkSmartPointer<vtkQuad>::New();

    for(uint i=0; i<connec.size(); i++)
    {
      for(uint k=0; k<connec[i].size(); k++)
        quad->GetPointIds()->SetId(k,connec[i][k]);

      cellArray->InsertNextCell(quad);
    }
    unstructuredGrid->SetCells(VTK_QUAD, cellArray);
  }
  else if (n_dim == 3 && n_node == 8)
  {
    vtkSmartPointer<vtkHexahedron> hex = 
      vtkSmartPointer<vtkHexahedron>::New();

    for(uint i=0; i<connec.size(); i++)
    {
      for(uint k=0; k<connec[i].size(); k++)
        hex->GetPointIds()->SetId(k,connec[i][k]);

      cellArray->InsertNextCell(hex);
    }
    unstructuredGrid->SetCells(VTK_HEXAHEDRON, cellArray);
  }
  else if (n_dim == 3 && n_node == 4)
  {
    vtkSmartPointer<vtkTetra> tet = 
      vtkSmartPointer<vtkTetra>::New();

    for(uint i=0; i<connec.size(); i++)
    {
      for(uint k=0; k<connec[i].size(); k++)
        tet->GetPointIds()->SetId(k,connec[i][k]);

      cellArray->InsertNextCell(tet);
    }
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);
  }
  else
  {
    throw std::runtime_error("Invalid element type in VtkIO.");
  }
}

void WriterVTK::write_material_id()
{
  vtkSmartPointer<vtkIntArray> material =
      vtkSmartPointer<vtkIntArray>::New();

  material->SetName("material");

  for(int i=0; i<mesh.get_n_elements(); i++)
    material->InsertNextValue( mesh.get_element_index(i) );

  unstructuredGrid->GetCellData()->AddArray(material);
}

void WriterVTK::write_point_data(const arma::vec & u, const std::string & name)
{
  const int usize = u.size();
  const double *sol = u.memptr();

  vtkSmartPointer<vtkDoubleArray> values =
      vtkSmartPointer<vtkDoubleArray>::New();

  values->SetName(name.c_str());

  for(int i=0; i<usize; i++)
    values->InsertNextValue( sol[i] );

  unstructuredGrid->GetPointData()->AddArray(values);
}

void WriterVTK::write_cell_data(const arma::vec & u, const std::string & name)
{
  const int usize = u.size();
  const double *sol = u.memptr();

  vtkSmartPointer<vtkDoubleArray> values =
      vtkSmartPointer<vtkDoubleArray>::New();

  values->SetName(name.c_str());

  for(int i=0; i<usize; i++)
    values->InsertNextValue( sol[i] );

  unstructuredGrid->GetCellData()->AddArray(values);
}

void WriterVTK::write_cell_vector(const std::vector<arma::vec3> & v, 
			      const std::string & name)
{
  const int vsize = v.size();

  vtkSmartPointer<vtkDoubleArray> vec =
      vtkSmartPointer<vtkDoubleArray>::New();
  vec->SetName(name.c_str());
  vec->SetNumberOfComponents(3);

  for(int i=0; i<vsize; i++)
  {
    double vd[3];
    vd[0] = v[i](0);
    vd[1] = v[i](1);
    vd[2] = v[i](2);
    vec->InsertNextTuple( vd );
  }
  unstructuredGrid->GetCellData()->AddArray(vec);
}

void WriterVTK::write_def_mesh(const std::vector<arma::vec3> & x, const string & file)
{
  set_points(x);
  set_cells();
  write_material_id();

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer->SetFileName(file.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();
  
  wfiles.push_back( file.c_str() );
}

void WriterVTK::write_def_mesh(const std::vector<arma::vec> & vu, 
			   const arma::vec & vms, const string & file)
{

  int nd = mesh.get_ndim();
  std::vector<arma::vec3> pts = mesh.get_points();

  // Compute deformed mesh
  for(uint i=0; i<pts.size(); i++)
    for(int d=0; d<nd; d++)
      pts[i](d) = pts[i](d) + vu[d][i];

  set_points(pts);
  set_cells();

  // Write cell data
  write_material_id();
  write_cell_data(vms, "vonmises");

  // Write displacements
  for(int d=0; d<nd; d++)
  {
    char aux[8];
    sprintf(aux,"u%d",(d+1));
    std::string dname(aux);
    write_point_data(vu[d],dname);
  }

  // Write file
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file.c_str() );
  writer->SetInputData(unstructuredGrid);
  writer->Write();
}

void WriterVTK::write_mesh(const string & file)
{
  std::vector<arma::vec3> pts = mesh.get_points();

  set_points(pts);
  set_cells();
  write_material_id();

  // Write file
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();
}

void WriterVTK::write_vtu(const std::string & file)
{
  // Write file
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();
}

void WriterVTK::write(const arma::vec & u, const string & file, bool binary)
{
  const std::vector<arma::vec3> & pts = mesh.get_points();
  
  set_points(pts);
  set_cells();
  write_point_data(u, "magnitude");
  write_material_id();

  // Write file
  if(binary)
  {
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
  }
  else // ascii
  {
    vtkSmartPointer<vtkUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(file.c_str());
    writer->SetFileTypeToASCII();
    writer->SetInputData(unstructuredGrid);
    writer->Write();
  }
  
  wfiles.push_back( file.c_str() );
}

*/