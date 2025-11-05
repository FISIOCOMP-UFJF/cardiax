#ifndef READMESH_H
#define	READMESH_H

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
#include <stdexcept>
#include <string>
#include <cstring>

//! Read mesh and return smart pointer to unstructured grid
vtkSmartPointer<vtkUnstructuredGrid> readmesh(const std::string & basename, 
                                              const std::string & fibname,
                                              int & fibermodel);

#endif	

