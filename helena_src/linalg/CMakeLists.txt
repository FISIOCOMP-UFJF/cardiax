# CmakeLists.txt in src dir
# Adds a library called LinAlg from the source files

link_directories(${AMGX_LIBRARIES})

include_directories (${AMGX_INCLUDE_DIRS})

add_library(LinAlg 
  tensor.cpp
  pardiso_solver.cpp
  petsc_vector.cpp 
  petsc_matrix.cpp 
  petsc_linear_solver.cpp
  )


target_link_libraries(LinAlg amgxsh gomp)
