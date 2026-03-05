set(CARDIAC_SRC
    cardiacproblem.cpp
    stimulus.cpp
    bidomain.cpp
    bidomain_bath.cpp
    bidomain_deformation.cpp
    monodomain.cpp
    monodomain_purkinje.cpp
    monodomain_deformation.cpp
    eikonal.cpp
    cardiac_cycle.cpp
    coupled_electromechanic.cpp
    coupled_purkinje.cpp
)
    
add_library(Cardiac ${CARDIAC_SRC})

target_include_directories(Cardiac PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Cardiac PUBLIC
    Util
    FEM
    Mesh
    ODE
    NLS
    LinAlg
    PDE
    ext_libs
)

