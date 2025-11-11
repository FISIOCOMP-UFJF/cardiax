# Lista de fontes (Isto está perfeito)
set(CARDIAC_SRC
    cardiacproblem.cpp
    stimulus.cpp
    bidomain.cpp
    bidomain_bath.cpp
    bidomain_deformation.cpp
    monodomain.cpp
    monodomain_purkinje.cpp
    monodomain_deformation.cpp
    cardiac_cycle.cpp
    #coupled_electromechanic.cpp
    coupled_purkinje.cpp
)
    
# Definição da biblioteca (Isto está perfeito)
add_library(Cardiac ${CARDIAC_SRC})

# [NOVO] Adiciona o diretório atual aos includes públicos.
# Isso permite que outros módulos (como 'app') façam #include "cardiacproblem.h"
target_include_directories(Cardiac PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# [MODIFICADO] Linka contra as bibliotecas
# Agora, linkamos contra os *alvos* internos do Cardiax
# e contra o alvo 'ext_libs' para TODAS as dependências externas.
target_link_libraries(Cardiac PUBLIC
    # Dependências Internas
    Util
    FEM
    Mesh
    ODE
    NLS
    LinAlg
    PDE

    # Dependências Externas (agrupadas)
    ext_libs
)

