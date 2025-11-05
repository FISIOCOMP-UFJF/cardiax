# Simple CMakeLists for building the applications
# Important to link against BLAS and LAPACK for good performance in Armadillo
# Bernardo M. Rocha

# Settings
# If using PARDISO, need to include pardiso500-GNU481-X86-64

set(CORELIBS blas hdf5 Util LinAlg ODE NLS Mesh FEM Cardiac PDE mpich)
set(BOOSTLIBS boost_system boost_filesystem)
set(PROGRAMS l2projection poisson elasticity nonlinearelas monodomain bidomain electromech)

# Build them all
foreach(program ${PROGRAMS})
  link_directories(${Boost_LIBRARIES})
  add_executable(${program} ${program}_solver.cpp)
  target_link_libraries(${program} ${CORELIBS} ${Boost_LIBRARIES} ${BOOSTLIBS} ${PETSC_LIBRARIES})
endforeach(program)

# Old way (para cada programa)
#add_executable(poisson poisson_solver.cpp)
#target_link_libraries(poisson ${CORELIBS} ${PETSC_LIBRARIES})
