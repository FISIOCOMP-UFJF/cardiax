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
  kerkoff2003.cpp
  ${agos_source} # Esta variável vem do seu CMakeLists.txt principal
)

target_include_directories(ODE PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(ODE PUBLIC
    Util
    LinAlg
    ext_libs # Nosso alvo com Armadillo, HDF5, PETSc, etc.
)

if(BUILD_TESTS)
    add_executable(testCellModel test/test_cellmodel.cpp)
    target_link_libraries(testCellModel PRIVATE ODE)
endif(BUILD_TESTS)