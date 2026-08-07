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

# [NOVO] Exporta os headers do LinAlg (tensor.h, etc.)
# para outros módulos (como FEM) que dependem dele.
target_include_directories(LinAlg PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# [MODIFICADO] Linka as dependências
# O alvo 'LinAlg' precisa de todas as dependências externas
# (PETSc, AMGX, Pardiso, OpenMP/gomp).
# Nós simplesmente linkamos contra o alvo 'ext_libs'
# que já contém TODAS elas.
target_link_libraries(LinAlg PUBLIC
    ext_libs
)