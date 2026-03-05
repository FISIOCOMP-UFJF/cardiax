# CmakeLists.txt in src dir
# Adds a library called nls from the source files

set (CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lz")

###link_directories ( ${ARMADILLO_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} )

link_directories(${AMGX_LIBRARIES})

include_directories (${AMGX_INCLUDE_DIRS})

add_library(NLS
    nonlinear_solver.cpp
    bfgs.cpp
    newton.cpp
    newton_ls.cpp
    test/simple.cpp
    test/test_simple.cpp
)

target_link_libraries(NLS amgxsh gomp)

#add_executable(testSimple test/test_simple.cpp test/simple.cpp newton.cpp nonlinear_solver.cpp)
#target_link_libraries(testSimple ${BLASVERSION} armadillo Util LinAlg ODE mpich ${PETSC_LIBRARIES})
