# CmakeLists.txt in src/ode dir
# Adds a library called ODE from the source files

add_library(ODE
  timestepper.cpp 
  cellmodel.cpp 
  cells.cpp
  ode_solver.cpp
  explicit_euler.cpp 
  runge_kutta4.cpp
  implicit_ode_solver.cpp
  implicit_euler.cpp
  simple_ode.cpp 
  fitz_hugh_nagumo.cpp
  nash_panfilov.cpp 
  mitchell_schaeffer.cpp
  mv.cpp
  ord.cpp
  luo_rudy.cpp 
  ten_tusscher2004.cpp
  ten_tusscher2006.cpp
  ten_tusscher_ta.cpp
  rice_ord.cpp
  rice_ten_tusscher.cpp
  mynash.cpp
  minimalmodel.cpp
  mm_silva.cpp
  torord_land.cpp
  ${agos_source} # Esta variável vem do seu CMakeLists.txt principal
)

# [NOVO] Exporta os headers do ODE
target_include_directories(ODE PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# [MODIFICADO] Linka as dependências de ODE
# Adicionamos 'ext_libs' porque os seus testes (testCellModel)
# precisam de Armadillo e HDF5, que estão em 'ext_libs'.
# Isso indica que a própria biblioteca 'ODE' (ou seus testes)
# usa bibliotecas externas.
target_link_libraries(ODE PUBLIC
    Util
    LinAlg
    ext_libs # Nosso alvo com Armadillo, HDF5, PETSc, etc.
)

# --- Programas de Teste ---
if(BUILD_TESTS)
    add_executable(testCellModel test/test_cellmodel.cpp)
    
    # [MODIFICADO] Linkagem do teste
    # Como 'ODE' já linka 'Util', 'LinAlg' e 'ext_libs'
    # como PUBLIC, o CMake cuida da linkagem transitiva.
    # 'testCellModel' recebe tudo o que precisa automaticamente.
    target_link_libraries(testCellModel PRIVATE ODE)

    # --- REMOVIDO ---
    # A linha antiga era:
    # target_link_libraries(testCellModel Util LinAlg ODE ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES} )
endif(BUILD_TESTS)