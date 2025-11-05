#include "helper.h"

vtkSmartPointer<vtkUnstructuredGrid> readmesh(const std::string & basename, 
                                              const std::string & fibname,
                                              int & fibermodel)
{
  std::string filename = basename;
  std::string fibrfile = fibname;
  std::string ptsname  = filename + ".pts";
  std::string elename  = filename + ".elem";

  // Create Unstructured Grid
  cout << "Creating unstructured grid" << endl;
  vtkSmartPointer<vtkUnstructuredGrid> ug =
          vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Read points
  int numPts;
  double px, py, pz;
  vtkSmartPointer<vtkPoints> points =
          vtkSmartPointer<vtkPoints>::New();

  ifstream ipts(ptsname.c_str());
  ipts >> numPts;
  for (int i = 0; i < numPts; i++)
  {
    ipts >> px >> py >> pz;
    points->InsertNextPoint(px, py, pz);
  }
  ipts.close();
  ug->SetPoints(points);

  // Read elements
  std::vector<int> ptnums;
  vtkSmartPointer<vtkCellArray> cellArray =
          vtkSmartPointer<vtkCellArray>::New();

  bool ishex = false;

  ifstream iele(elename.c_str());
  int numElem;
  std::string ecode;
  int n1, n2, n3, n4, n5, n6, n7, n8, temp;
  iele >> numElem;
  for (int i = 0; i < numElem; i++)
  {
    iele >> ecode;

    if (ecode == "Hx")
    {
      ishex = true;
      vtkSmartPointer<vtkHexahedron> hex = 
              vtkSmartPointer<vtkHexahedron>::New();
      
      iele >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> n8 >> temp;
      
      hex->GetPointIds()->SetId(0, n1);
      hex->GetPointIds()->SetId(1, n2);
      hex->GetPointIds()->SetId(2, n3);
      hex->GetPointIds()->SetId(3, n4);
      hex->GetPointIds()->SetId(4, n5);
      hex->GetPointIds()->SetId(5, n6);
      hex->GetPointIds()->SetId(6, n7);
      hex->GetPointIds()->SetId(7, n8);
      cellArray->InsertNextCell(hex);
    }    
    else
      throw std::runtime_error("Element type not available.");    
  }

  if (ecode == "Hx")
  {
    ug->SetCells(VTK_HEXAHEDRON, cellArray);
    cout << "Element type: HEX" << endl;
  }
  iele.close();

  // Print information
  int ncells = ug->GetNumberOfCells();
  int npoints = ug->GetNumberOfPoints();
  cout << "Number of cells: " << ncells << endl;
  cout << "Number of points: " << npoints << endl;

  // Read fiber direction information
  cout << "Reading fiber information" << endl;
  ifstream ifib(fibrfile.c_str());

  vtkSmartPointer<vtkDoubleArray> f0 =
          vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> s0 =
          vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> n0 =
          vtkSmartPointer<vtkDoubleArray>::New();

  f0->SetName("Fiber");
  s0->SetName("Sheet");
  n0->SetName("Normal");
  f0->SetNumberOfComponents(3);
  s0->SetNumberOfComponents(3);
  n0->SetNumberOfComponents(3);

  int fibtype;
  ifib >> fibtype;

  cout << "Fiber model: " << fibtype << " (0-ortotropic, 1-fiber only)" << endl;
  for (int i = 0; i < npoints; i++)
  {
    double vecf[3] = {0}, vecs[3] = {0}, vecn[3] = {0};

    ifib >> vecf[0] >> vecf[1] >> vecf[2]; // fiber
    if(fibtype == 0)
    {
        ifib >> vecs[0] >> vecs[1] >> vecs[2]; // sheet
        ifib >> vecn[0] >> vecn[1] >> vecn[2]; // normal
    }
    f0->InsertNextTuple(vecf);
    s0->InsertNextTuple(vecs);
    n0->InsertNextTuple(vecn);
  }
  ifib.close();

  ug->GetPointData()->AddArray(f0);
  ug->GetPointData()->AddArray(s0);
  ug->GetPointData()->AddArray(n0);   
  
  fibermodel = fibtype;

  return ug;
}
