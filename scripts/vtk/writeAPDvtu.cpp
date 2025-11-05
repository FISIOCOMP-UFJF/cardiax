#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <iostream>
#include <fstream>
 
using namespace std;

int main (int argc, char *argv[])
{
  // parse command line arguments
  if(argc < 4)
    {
      cout << endl << " Usage: " << argv[0] << " <filename>.vtu";
      cout << " celltypes>.typ <apdfile>.dat [actfile.dat] [repfile.dat]" ;
      cout << endl << endl;
      return EXIT_FAILURE;
    }
 
  std::string filename, ctypfile, apdmfile, actfile, repfile;
  filename = argv[1];  
  ctypfile = argv[2];
  apdmfile = argv[3];  
  if(argc==6)
  {
    actfile = argv[4];
    repfile = argv[5];
  }   
 
  // read all the data from the file
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkSmartPointer<vtkUnstructuredGrid> ug =
     vtkSmartPointer<vtkUnstructuredGrid>::New();
  ug->ShallowCopy(reader->GetOutput());
  
  // information
  int ncells  = ug->GetNumberOfCells();
  int npoints = ug->GetNumberOfPoints();
  cout << "Number of cells: " << ncells << endl;
  cout << "Number of points: " << npoints << endl;
  
  //
  // read cell phenotypes
  //
  cout << "Reading cell types" << endl;
  int num;
  ifstream inct;
  inct.open(ctypfile.c_str());
  inct >> num;

  vtkSmartPointer<vtkIntArray> types =
    vtkSmartPointer<vtkIntArray>::New();
  types->SetNumberOfComponents(1);
  types->SetName("CellType");

  for(int i=0; i<num; i++)
  {
    int temp;
	inct >> temp;
	types->InsertNextValue(temp);
  }
  inct.close();

  //
  // append cell phenotypes to unstructured grid
  //
  cout << "Writing cell types to meshAPDMapColor.vtu" << endl;
  ug->GetPointData()->AddArray(types);

  //
  // read APD map color
  //
  cout << "Reading APD90 data" << endl;
  vtkSmartPointer<vtkFloatArray> values =
    vtkSmartPointer<vtkFloatArray>::New();
  values->SetNumberOfComponents(1);
  values->SetName("APD");

  ifstream inapd(apdmfile.c_str());
  for(int i=0; i<npoints; i++)
  {
    float temp;
    inapd >> temp;
    values->InsertNextValue(temp);
  }
  inapd.close();

  // append APD data to Unstructured Grid
  cout << "Writing APD90 data" << endl;
  ug->GetPointData()->AddArray(values);
  
  if(argc==6)
  {
    cout << "Reading ACTivation data" << endl;
    vtkSmartPointer<vtkFloatArray> valuesACT =
      vtkSmartPointer<vtkFloatArray>::New();
    valuesACT->SetNumberOfComponents(1);
    valuesACT->SetName("ACT");

    ifstream inact(actfile.c_str());
    for(int i=0; i<npoints; i++)
    {
      float temp;
      inact >> temp;
      valuesACT->InsertNextValue(temp);
    }
    inact.close();   
    cout << "Writing ACT data" << endl;
    ug->GetPointData()->AddArray(valuesACT);
    
    cout << "Reading REPolarization data" << endl;
    vtkSmartPointer<vtkFloatArray> valuesREP =
      vtkSmartPointer<vtkFloatArray>::New();
    valuesREP->SetNumberOfComponents(1);
    valuesREP->SetName("REP");

    ifstream inrep(repfile.c_str());
    for(int i=0; i<npoints; i++)
    {
      float temp;
      inrep >> temp;
      valuesREP->InsertNextValue(temp);
    }
    inrep.close();   
    cout << "Writing REP data" << endl;
    ug->GetPointData()->AddArray(valuesREP);
  }
  
  //
  // write VTU + APD Map color
  //
  cout << "Writing meshAPDMapColor.vtu" << endl;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName("meshAPDMapColor.vtu");
  writer->SetInput(ug);
  writer->Write();
  
  cout << "Done" << endl;
 
  return EXIT_SUCCESS;
}
