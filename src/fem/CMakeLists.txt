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

# [NOVO] Exporta os headers do FEM para outros alvos
# que linkarem contra ele (como o 'cardiax_app')
target_include_directories(FEM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
  
# [MODIFICADO] Linka as dependências de FEM
# PUBLIC significa que qualquer alvo que linkar 'FEM'
# também receberá 'Util', 'LinAlg' e 'ext_libs'.
target_link_libraries(FEM PUBLIC
    Util
    LinAlg
    ext_libs # Adiciona nosso alvo de interface
)

# --- Testes ---

# [SEM MUDANÇA] Define o executável de teste
add_executable(testQuadrature test/test_quadrature.cpp)

# [MODIFICADO] Linkagem do teste
# Esta é a principal vantagem do CMake moderno.
# Como 'testQuadrature' depende de 'FEM', e 'FEM'
# já linka 'Util', 'LinAlg' e 'ext_libs' como PUBLIC,
# não precisamos repetir nada. O CMake gerencia a
# linkagem transitiva.
target_link_libraries(testQuadrature PRIVATE FEM)

# --- REMOVIDO ---
# A linha antiga foi:
# target_link_libraries(testQuadrature FEM Util LinAlg ${PETSC_LIBRARIES})