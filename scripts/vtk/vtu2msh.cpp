
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
 
int main (int argc, char *argv[])
{
  //parse command line arguments
  if(argc != 2)
    {
    std::cerr << "Usage: " << argv[0]
              << " filename.vtu" << std::endl;
    return EXIT_FAILURE;
    }
 
  std::string filename = argv[1];
 
  //read all the data from the file
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkSmartPointer<vtkUnstructuredGrid> ug =
     vtkSmartPointer<vtkUnstructuredGrid>::New();
  ug->ShallowCopy(reader->GetOutput());
  
  // imprime informacoes
  int ncells  = ug->GetNumberOfCells();
  int npoints = ug->GetNumberOfPoints();

  cout << "Number of cells: " << ncells << endl;
  cout << "Number of points: " << npoints << endl;

  // escreve arquivo .msh
  ofstream out("ftest.msh");

  out << "$MeshFormat" << endl;
  out << "2.0 0 8" << endl;
  out << "$EndMeshFormat" << endl;

  vtkSmartPointer<vtkPoints> pts = ug->GetPoints();

  out << "$Nodes" << endl;
  out << npoints << endl;
  for(int i=0; i<pts->GetNumberOfPoints(); i++)
  {
    double *pt = pts->GetPoint(i);
    out << i+1 << " " << pt[0] << " " << pt[1] << " " << pt[2] << endl;
  }
  out << "$EndNodes" << endl;

  out << "$Elements" << endl;
  out << ncells << endl;
  for(int i=0; i<ncells; i++)
  {    
    vtkSmartPointer<vtkCell> cell = ug->GetCell(i);
    vtkSmartPointer<vtkIdList> pids = cell->GetPointIds();

    out << i+1 << " 5 3 1 10 99 ";
    for(int j=0; j<pids->GetNumberOfIds(); j++)
    {
      out << cell->GetPointId(j) + 1 << " ";
    }
    out << endl;
  }
  out << "$EndElements" << endl;
  out.close();
 
  return EXIT_SUCCESS;
}
