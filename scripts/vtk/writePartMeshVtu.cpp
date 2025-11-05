/*
 * writePartMeshVtu.cpp
 *
 *  Created on: 29/10/2013
 *      Author: rocha
 */

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
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
      std::cout << "Usage: <basename> <file.part>" << endl;
      return EXIT_FAILURE;
    }

  std::string filename = argv[1];
  std::string partfile = argv[2];

  std::string ptsname = filename + ".pts";
  std::string elename = filename + ".elem";

  // Create Unstructured Grid
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
	//cout << px << "\t" << py << "\t" << pz << endl;
  }
  ipts.close();
  ug->SetPoints(points);

  // Read elements
  std::vector<int> ptnums;

  vtkSmartPointer<vtkCellArray> cellArray =
		  vtkSmartPointer<vtkCellArray>::New();

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

  }
  ug->SetCells(VTK_HEXAHEDRON, cellArray);
  iele.close();

  // Read partitions
  ifstream ipart(partfile.c_str());

  vtkSmartPointer<vtkIntArray> material =
      vtkSmartPointer<vtkIntArray>::New();

  material->SetName("Subdomain");
  int id;
  for(int i=0; i<numElem; i++)
  {
	ipart >> id;
	//cout << id << endl;
    material->InsertNextValue(id);
  }
  ipart.close();
  ug->GetCellData()->AddArray(material);

  // imprime informacoes
  int ncells  = ug->GetNumberOfCells();
  int npoints = ug->GetNumberOfPoints();
  cout << "Number of cells: " << ncells << endl;
  cout << "Number of points: " << npoints << endl;

  //
  // Write VTU + Partition
  //
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName("PartitionedMesh.vtu");
  writer->SetInput(ug);
  writer->Write();

  cout << "Done" << endl;
  return EXIT_SUCCESS;
}



