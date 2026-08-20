# CmakeLists.txt in src/fem dir
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

target_include_directories(FEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(FEM PUBLIC
    Util
    LinAlg
    ext_libs 
)

add_executable(testQuadrature test/test_quadrature.cpp)

target_link_libraries(testQuadrature PUBLIC FEM)
