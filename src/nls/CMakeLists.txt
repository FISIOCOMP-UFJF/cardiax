# --- REMOVIDO ---
# set (CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lz")
# link_directories(...)
# include_directories(...)
# Todas as linhas acima são desnecessárias. O alvo 'ext_libs'
# cuida de AMGX, zlib, PETSc, etc.

# Define a biblioteca NLS (sem os arquivos de teste)
add_library(NLS
    nonlinear_solver.cpp
    bfgs.cpp
    newton.cpp
    newton_ls.cpp
)

# [NOVO] Exporta os headers do NLS
target_include_directories(NLS PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# [MODIFICADO] Linka as dependências de NLS
# 'NLS' (Non-Linear Solver) provavelmente depende de 'LinAlg'
# e das bibliotecas externas.
target_link_libraries(NLS PUBLIC
    LinAlg
    ext_libs # Fornece PETSc, AMGX, gomp, zlib, etc.
)

# --- Teste ---
# (Baseado nas suas linhas comentadas)

# [MODIFICADO] Define o executável de teste
# Ele compila os fontes de teste
add_executable(testSimple EXCLUDE_FROM_ALL
    test/simple.cpp
    test/test_simple.cpp
)

# [MODIFICADO] Linka o teste
# O teste precisa da biblioteca NLS que ele está testando.
# Ele também parecia precisar de 'ODE' e 'Util',
# baseado no seu comentário original.
target_link_libraries(testSimple PRIVATE
    NLS
    ODE
    Util
)
# Nota: Como NLS linka 'LinAlg' e 'ext_libs' como PUBLIC,
# 'testSimple' recebe todas essas dependências transitivamente.
# Não precisamos repetir ${PETSC_LIBRARIES}, armadillo, etc.