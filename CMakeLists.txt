# CARDIAX
# author: Bernardo Martins Rocha
# mail: bernardomartinsrocha@gmail.com

cmake_minimum_required (VERSION 2.8)
project (Cardiax)

# Configure PETSc dir and arch
#set (ENV{PETSC_DIR} "$ENV{PETSC_DIR}")
#set (ENV{PETSC_ARCH} "$ENV{PETSC_ARCH}")
#find_path (PETSC_DIR include/petsc.h
#  HINTS ENV PETSC_DIR
#  PATHS
#  ${CMAKE_SYSTEM_PREFIX_PATH}
#  DOC "PETSc Directory")
#set(PETSC_ARCH 'arch-linux2-c-debug')

# Configure cmake modules dir
set (CMAKE_C_FLAGS "-O3 -Wall -Wno-deprecated")
set (CMAKE_CXX_FLAGS "-O3 -Wall -Wno-deprecated")
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

set (MPI_C_LIBRARIES "/home/joventino/source/petsc-3.19.6/arch-linux2-c-debug/lib")
set (MPI_C_INCLUDE_PATH "/home/joventino/source/petsc-3.19.6/arch-linux2-c-debug/include")
set (MPI_CXX_LIBRARIES "/home/joventino/source/petsc-3.19.6/arch-linux2-c-debug/lib")
set (MPI_CXX_INCLUDE_PATH "/home/joventino/source/petsc-3.19.6/arch-linux2-c-debug/include")
set (MPI_C_COMPILER "/home/joventino/source/petsc-3.19.6/arch-linux2-c-debug/bin/mpicc")
set (MPI_CXX_COMPILER "/home/joventino/source/petsc-3.19.6/arch-linux2-c-debug/bin/mpicxx")

set (MPI_lib_LIBRARY "/home/joventino/source/petsc-3.19.6/arch-linux2-c-debug/lib")

#set (HDF5_INCLUDE_DIRS "/opt/ohpc/pub/libs/gnu12/mpich/hdf5/1.10.8/include/")
#set (HDF5_LIBRARIES "/opt/ohpc/pub/libs/gnu12/mpich/hdf5/1.10.8/lib/")

set (HDF5_INCLUDE_DIRS "/home/joventino/source/hdf5/include/")
set (HDF5_LIBRARIES "/home/joventino/source/hdf5/lib/")

set (Boost_INCLUDE_DIRS "/home/joventino/source/boost_1_81_0/boost/")
set (Boost_LIBRARIES "/home/joventino/source/boost_1_81_0/stage/lib")

#set (BLAS_INCLUDE_DIRS "/usr/include/")
set (BLAS_LIBRARIES "/usr/lib64/mpich/lib/")

set (ARMADILLO_LIBRARY "/home/joventino/source/armadillo/usr/lib64/")
set (ARMADILLO_INCLUDE_DIR "/home/joventino/source/armadillo/usr/include/")

set (AMGX_INCLUDE_DIRS "/home/joventino/source/amgx/include/")
set (AMGX_LIBRARIES "/home/joventino/source/amgx/lib/")
set (AMGX_DIR "/home/joventino/source/amgx")

# Find external packages
find_package(MPI REQUIRED)
find_package(PETSc REQUIRED)
#find_package(HDF5 REQUIRED)
#find_package(Boost REQUIRED)
find_package(Armadillo REQUIRED)
#find_package(AMGX REQUIRED)

message(${ARMADILLO_LIBRARY})
message(${HDF5_LIBRARIES})
message(${HDF5_LIBRARIES})
#message(${Boost_LIBRARIES})
message(${AMGX_LIBRARIES})

# Includes and libraries for external packages
#set (PARDISO_LIB /usr/lib/pardiso)


# OLD VTK stuff
#include (${CMAKE_ROOT}/Modules/FindVTK.cmake)
#include (${PROJECT_SOURCE_DIR}/cmake/Modules/FindVTK.cmake)
###if (USE_VTK_FILE)
###   include(${USE_VTK_FILE})
###endif (USE_VTK_FILE)
###find_package(VTK COMPONENTS REQUIRED)

link_directories( 
  ${MPI_CXX_LIBRARIES}
  ${PETSC_LIBRARIES} 
  ${ARMADILLO_LIBRARIES} 
  ${HDF5_LIBRARIES}
  ${Boost_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${AMGX_LIBRARIES}
  )

include_directories (
  ${PROJECT_SOURCE_DIR}/src/  
  ${PROJECT_SOURCE_DIR}/src/util/ 
  ${PROJECT_SOURCE_DIR}/src/linalg/ 
  ${PROJECT_SOURCE_DIR}/src/odes/
  ${PROJECT_SOURCE_DIR}/src/fem/
  ${PROJECT_SOURCE_DIR}/src/pde/
  ${Boost_INCLUDE_DIRS}
  ${PETSC_INCLUDES}
  ${HDF5_INCLUDE_DIRS}
  ${ARMADILLO_INCLUDE_DIRS}
  ${AMGX_INCLUDE_DIRS} 
  )

# Configure some options for compiling or not extra packages, etc... 
option(BUILD_TESTS "Build tests" ON)     
option(USE_AGOS "Use AGOS for ODE solvers" OFF)
  
add_definitions(${PETSC_DEFINITIONS})
add_subdirectory(src/util)
add_subdirectory(src/linalg)
add_subdirectory(src/mesh)
add_subdirectory(src/fem)
add_subdirectory(src/odes)
add_subdirectory(src/pdes)
add_subdirectory(src/cardiac)
add_subdirectory(src/nls)
add_subdirectory(app)

if(BUILD_TESTS)
  add_subdirectory(test)   
endif(BUILD_TESTS)

#
#message(STATUS "OK Listing all include directories ")
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS " dir=${dir}")
#endforeach()  

#message(STATUS "OK Listing all library directories ")
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS " dir=${dir}")
#endforeach()  

