# CmakeLists.txt in src/pde dir
# Adds a library called PDE from the source files

set(
  PDE_SRC 
  poisson.cpp 
  laplace.cpp 
  ###mortarfem.cpp
  elasticity.cpp 
  mech_utils.cpp
  nonlinear_elasticity.cpp
  total_lagrangian.cpp
  #total_lagrangian_snes.cpp
  updated_lagrangian.cpp 
  materials/hyperelastic_material.cpp
  materials/incompressible_material.cpp
  materials/material_data.cpp
  materials/saint_venant_kirchoff.cpp
  materials/compressible_neo_hookean.cpp
  materials/neo_hookean.cpp
  materials/guccione.cpp
  materials/exponential.cpp
  materials/holzapfel_ogden.cpp
)

# Definição da biblioteca (Isto está perfeito)
add_library(PDE ${PDE_SRC})

# [NOVO] Exporta os headers do PDE para outros módulos
target_include_directories(PDE PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# [MODIFICADO] Linka as dependências como PUBLIC
# Isso garante que qualquer alvo que linkar 'PDE' (como 'Cardiac')
# também receba as propriedades de 'FEM', 'ODE', e 'LinAlg'.
target_link_libraries(PDE PUBLIC
    FEM
    ODE
    LinAlg
)