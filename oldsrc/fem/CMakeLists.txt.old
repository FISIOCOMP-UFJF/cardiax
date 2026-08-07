# CmakeLists.txt in src dir
# Adds a library called FEM from the source files

add_library(FEM 
  fe.cpp
  mixed_fe.cpp
  fespace.cpp
  mixed_fespace.cpp
  hermite_fespace.cpp
  fedata.cpp
  fetools.cpp
  #../mesh/gmshio.cpp
  quadrature.cpp
  mapping.cpp
  #../mesh/mesh.cpp
  #../mesh/box_mesh.cpp
  #../mesh/purkinje_mesh.cpp
  projection.cpp
  #../mesh/writer.cpp
  #../mesh/writer_vtk.cpp
  #../mesh/writer_hdf5.cpp
  )
  
target_link_libraries(FEM Util LinAlg armadillo hdf5)

#
# Binaries for tests
#
add_executable(testQuadrature test/test_quadrature.cpp)
target_link_libraries(testQuadrature FEM Util LinAlg armadillo ${PETSC_LIBRARIES})