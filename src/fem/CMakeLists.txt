# CmakeLists.txt in src dir
# Adds a library called FEM from the source files

add_library(FEM 
  fe.cpp
  fespace.cpp
  fedata.cpp
  fetools.cpp  
  mixed_fe.cpp
  mixed_fespace.cpp
  hermite_fespace.cpp
  quadrature.cpp
  mapping.cpp
  projection.cpp  
  )
  
target_link_libraries(FEM Util LinAlg)

#
# Binaries for tests
#
add_executable(testQuadrature test/test_quadrature.cpp)
target_link_libraries(testQuadrature FEM Util LinAlg ${PETSC_LIBRARIES})
