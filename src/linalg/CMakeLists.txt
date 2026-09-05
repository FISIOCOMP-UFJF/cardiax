# --- REMOVIDO ---
# link_directories(${AMGX_LIBRARIES})
# include_directories (${AMGX_INCLUDE_DIRS})
# Estas diretivas são globais e devem ser evitadas.
# As informações de include e link do AMGX já estão no
# alvo 'ext_libs' que definimos no CMakeLists.txt principal.

# Definição da biblioteca (Isto está perfeito)
add_library(LinAlg 
  tensor.cpp
  pardiso_solver.cpp
  petsc_vector.cpp 
  petsc_matrix.cpp 
  petsc_linear_solver.cpp
)

target_include_directories(LinAlg PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    $ENV{CONDA_PREFIX}/include
    /usr/local/cuda/include
)


target_link_libraries(LinAlg PUBLIC
    ext_libs
)