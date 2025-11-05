#include <cassert>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

int main (int argc, char *argv[])
{
  // parse command line arguments
  if(argc != 3)
  {
    std::cerr << "\n Usage: " << argv[0]
              << " fenicsFileName.vtu myFilename.vtu\n" << std::endl;
    return EXIT_FAILURE;
  }
  
  int ncells1, ncells2, npoints1, npoints2;
  std::string file1 = argv[1];
  std::string file2 = argv[2];

  // read all the data from the file1
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader1 =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader1->SetFileName(file1.c_str());
  reader1->Update();

  vtkSmartPointer<vtkUnstructuredGrid> ug1 =
     vtkSmartPointer<vtkUnstructuredGrid>::New();
  ug1->ShallowCopy(reader1->GetOutput());
  
  // imprime informacoes
  ncells1  = ug1->GetNumberOfCells();
  npoints1 = ug1->GetNumberOfPoints();

  cout << "File 1" << endl;
  cout << "Number of cells: " << ncells1 << endl;
  cout << "Number of points: " << npoints1 << endl;

  // read all data from the file2
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader2 =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader2->SetFileName(file2.c_str());
  reader2->Update();

  vtkSmartPointer<vtkUnstructuredGrid> ug2 =
     vtkSmartPointer<vtkUnstructuredGrid>::New();
  ug2->ShallowCopy(reader2->GetOutput());
  
  // imprime informacoes
  ncells2  = ug2->GetNumberOfCells();
  npoints2 = ug2->GetNumberOfPoints();

  cout << "File 2" << endl;
  cout << "Number of cells: " << ncells2 << endl;
  cout << "Number of points: " << npoints2 << endl;

  assert(ncells1 == ncells2);
  assert(npoints1 == npoints2);
  
  cout << "Comparing displacements" << endl;

  vtkSmartPointer<vtkDoubleArray> array1 =
	  vtkDoubleArray::SafeDownCast(ug1->GetPointData()->GetVectors("f_2"));	    
   
  vtkSmartPointer<vtkDoubleArray> array2x =
	  vtkDoubleArray::SafeDownCast(ug2->GetPointData()->GetScalars("DisplacementX"));
  vtkSmartPointer<vtkDoubleArray> array2y =
	  vtkDoubleArray::SafeDownCast(ug2->GetPointData()->GetScalars("DisplacementY"));
  vtkSmartPointer<vtkDoubleArray> array2z =
	  vtkDoubleArray::SafeDownCast(ug2->GetPointData()->GetScalars("DisplacementZ"));
      
  vtkSmartPointer<vtkDoubleArray> arrayDifX = 
	  vtkSmartPointer<vtkDoubleArray>::New(); 

  vtkSmartPointer<vtkDoubleArray> arrayDifY = 
	  vtkSmartPointer<vtkDoubleArray>::New(); 

  vtkSmartPointer<vtkDoubleArray> arrayDifZ = 
	  vtkSmartPointer<vtkDoubleArray>::New(); 
  
  for(int i=0; i<npoints1; i++)
  {
    double* u1 = array1->GetTuple(i);  
    double u2x = array2x->GetValue(i);
    double u2y = array2y->GetValue(i);
    double u2z = array2z->GetValue(i);
    
    arrayDifX->InsertValue(i,fabs(u1[0]-u2x));
    arrayDifY->InsertValue(i,fabs(u1[1]-u2y));
    arrayDifZ->InsertValue(i,fabs(u1[2]-u2z));
    
    cout << std::scientific << u1[0] << "\t" << u1[1] << "\t" << u1[2] << "\t";
    cout << std::scientific << u2x << "\t" << u2y << "\t" << u2z;
    cout << endl;
        
    //cout << std::scientific << fabs(u1[0]-u2x) << "\t";
    //cout << fabs(u1[1]-u2y) << "\t" << fabs(u1[2]-u2z) << endl; 
	  
  }
  
  double range[2];
  
  cout << "Errors in displacement" << endl;
  
  arrayDifX->GetRange(range);
  cout << " max in X: " << range[1] << endl;
  
  arrayDifY->GetRange(range);  
  cout << " max in Y: " << range[1] << endl;
  
  arrayDifZ->GetRange(range);
  cout << " max in Z: " << range[1] << endl;

  return EXIT_SUCCESS;
}
