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
#include <vtkPointDataToCellData.h>
#include <iostream>
#include <fstream>

#include "helper.h"

using namespace std;

int main(int argc, char *argv[])
{
  //parse command line arguments
  if (argc != 3)
  {
    std::cout << "Usage: <basename> <file.fib>" << endl;
    return EXIT_FAILURE;
  }

  int fibtype;
  std::string filename = argv[1];
  std::string fibrfile = argv[2];
  
  // read mesh and create unstructured grid
  vtkSmartPointer<vtkUnstructuredGrid> ug = readmesh(filename, fibrfile, fibtype);
    
  //
  // Write fiber directions
  //
  cout << "Writing mesh: MeshFibers.vtu" << endl;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
          vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer->SetFileName("MeshFibers.vtu");
  writer->SetInput(ug);
  writer->Write();

  //
  // Convert Fibers from Points to Cells
  //
  vtkSmartPointer<vtkXMLUnstructuredGridReader> ugreader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  ugreader->SetFileName("MeshFibers.vtu");

  vtkSmartPointer<vtkPointDataToCellData> p2c = 
    vtkSmartPointer<vtkPointDataToCellData>::New();
  p2c->SetInputConnection(ugreader->GetOutputPort());

  cout << "Writing mesh: MeshFibersCellData.vtu" << endl;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> cwriter =
          vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  cwriter->SetFileName("MeshFibersCellData.vtu");
  cwriter->SetInputConnection(p2c->GetOutputPort());
  cwriter->Write();

  //
  // Read the cell data and write to plain text file 
  //
  vtkSmartPointer<vtkXMLUnstructuredGridReader> ugc_reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  ugc_reader->SetFileName("MeshFibersCellData.vtu");
  ugc_reader->Update();

  vtkSmartPointer<vtkUnstructuredGrid> ugCell =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  ugCell->ShallowCopy(ugc_reader->GetOutput());

  vtkSmartPointer<vtkDoubleArray> afibersCells =
    vtkDoubleArray::SafeDownCast(ugCell->GetCellData()->GetVectors("Fiber"));

  vtkSmartPointer<vtkDoubleArray> asheetsCells =
    vtkDoubleArray::SafeDownCast(ugCell->GetCellData()->GetVectors("Sheet"));

  vtkSmartPointer<vtkDoubleArray> anormalsCells =
    vtkDoubleArray::SafeDownCast(ugCell->GetCellData()->GetVectors("Normal"));
  
  cout << "Writing plain text file: MeshFibersCellData.fib" << endl;
  ofstream fibfile("MeshFibersCellData.fib");  
  
  if(fibtype==0) fibfile << "0" << endl;
  
  for(int i=0; i<afibersCells->GetNumberOfTuples(); i++)
  {
    double * uf = afibersCells->GetTuple(i);
    fibfile << uf[0] << " " << uf[1] << " " << uf[2] << " ";
    if(fibtype == 0)
    {
      double * us = asheetsCells->GetTuple(i);
      double * un = anormalsCells->GetTuple(i);
      fibfile << us[0] << " " << us[1] << " " << us[2] << " ";
      fibfile << un[0] << " " << un[1] << " " << un[2] << " ";
    }
    fibfile << endl;
  }
  fibfile.close();

  cout << "Done" << endl;

  return EXIT_SUCCESS;
}
