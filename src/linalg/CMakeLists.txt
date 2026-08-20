add_library(LinAlg 
  tensor.cpp
  pardiso_solver.cpp
  petsc_vector.cpp 
  petsc_matrix.cpp 
  petsc_linear_solver.cpp
)

target_include_directories(LinAlg PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(LinAlg PUBLIC
    ext_libs
)