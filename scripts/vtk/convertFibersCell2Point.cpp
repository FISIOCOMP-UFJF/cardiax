#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkHexahedron.h>
#include <iostream>
#include <fstream>

using namespace std;

int main (int argc, char *argv[])
{
  //parse command line arguments
  if(argc != 3)
  {
    std::cout << "Usage: <basename> <file.fib>" << endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];
  std::string fibrfile = argv[2];

  std::string ptsname = filename + ".pts";
  std::string elename = filename + ".elem";

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
  for(int i=0; i<numPts; i++)
  {
    ipts >> px >> py >> pz;
    points->InsertNextPoint(px,py,pz);
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
  int n1,n2,n3,n4,n5,n6,n7,n8,temp;
  iele >> numElem;
  for(int i=0; i<numElem; i++)
  {
    iele >> ecode;

    if(ecode == "Hx")
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
    else if(ecode == "Tt")
    {
      vtkSmartPointer<vtkTetra> tet =
    		  vtkSmartPointer<vtkTetra>::New();
      iele >> n1 >> n2 >> n3 >> n4 >> temp;
      tet->GetPointIds()->SetId(0, n1);
      tet->GetPointIds()->SetId(1, n2);
      tet->GetPointIds()->SetId(2, n3);
      tet->GetPointIds()->SetId(3, n4);
      cellArray->InsertNextCell(tet);
    }
  }

  if(ecode == "Hx")
    ug->SetCells(VTK_HEXAHEDRON, cellArray);
  else if(ecode == "Tt")
    ug->SetCells(VTK_TETRA, cellArray);  
  iele.close();

  // Print information
  int ncells  = ug->GetNumberOfCells();
  int npoints = ug->GetNumberOfPoints();
  cout << "Number of cells: " << ncells << endl;
  cout << "Number of points: " << npoints << endl;

  // Read fiber direction information
  cout << "Reading FIBER information" << endl;
  ifstream ifib(fibrfile.c_str());

  vtkSmartPointer<vtkDoubleArray> f0 = vtkSmartPointer<vtkDoubleArray>::New();

  f0->SetName("Fiber");
  f0->SetNumberOfComponents(3);

  int fibtype;
  ifib >> fibtype;

  for(int i=0; i<ncells; i++)
  {
    double vecf[3];
    ifib >> vecf[0] >> vecf[1] >> vecf[2]; // fiber    
    f0->InsertNextTuple(vecf);
  }
  ifib.close();

  ug->GetCellData()->AddArray(f0);

  // CONVERT FROM CELL TO POINT  

  //
  // Write fiber directions
  //
  cout << "Writing mesh with fiber, sheet and normal directions" << endl;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer->SetFileName("NewMeshWithFibers.vtu");
  writer->SetInput(ug);
  writer->Write();
  cout << "Done" << endl;

  return EXIT_SUCCESS;
}
