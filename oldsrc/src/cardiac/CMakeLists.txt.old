# CmakeLists.txt in src dir
# Adds a library called problems from the source files

link_directories ( ${ARMADILLO_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ) #${VTK_LIBRARY_DIRS}

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
	
add_library(Cardiac ${CARDIAC_SRC})
target_link_libraries(Cardiac Util FEM Mesh ODE NLS LinAlg armadillo PDE) #${VTK_LIBRARIES})
