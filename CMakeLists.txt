# CARDIAX
# author: Bernardo Martins Rocha
# mail: bernardomartinsrocha@gmail.com

cmake_minimum_required(VERSION 3.12) 
project(Cardiax C CXX)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 OLD)
endif()

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type: Debug, Release, RelWithDebInfo, MinSizeRel")
  message(STATUS "No build type specified, defaulting to 'Release'.")
  message(STATUS "Use -DCMAKE_BUILD_TYPE=Debug for debugging.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

if (NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "The environment variable PETSC_DIR must be set to the PETSc directory.")
endif()
if (NOT DEFINED ENV{PETSC_ARCH})
    message(WARNING "The environment variable PETSC_ARCH is not set. Using 'arch-linux2-c-debug' as the default.")
    set(ENV{PETSC_ARCH} "arch-linux2-c-debug")
endif()

find_package(PETSc REQUIRED)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
# find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ScaLAPACK REQUIRED)
find_library(ARMADILLO_LIBRARY_PATH armadillo HINTS $ENV{ARMADILLO_ROOT} PATH_SUFFIXES lib lib64)


set(AMGX_INCLUDE_DIRS $ENV{AMGX_ROOT}/include)
# set(AMGX_LIBRARIES $ENV{AMGX_ROOT}/lib/libamgx.a) 
find_library(AMGX_LIBRARY_PATH amgx HINTS $ENV{AMGX_ROOT} PATH_SUFFIXES lib lib64)

find_library(CUDART_LIB cudart HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
find_library(CUSPARSE_LIB cusparse HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
find_library(CUSOLVER_LIB cusolver HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
find_library(NVTX_LIB nvToolsExt HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
find_library(CUBLAS_LIB cublas HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)


set(AMGX_LIBRARIES ${AMGX_LIBRARY_PATH})
set(CUDA_LIBRARIES ${CUDART_LIB} ${CUSPARSE_LIB} ${CUSOLVER_LIB} ${NVTX_LIB} ${CUBLAS_LIB})


add_library(ext_libs INTERFACE)

target_include_directories(ext_libs INTERFACE
    ${PROJECT_SOURCE_DIR}/src
    ${MPI_C_INCLUDE_DIRS}
    ${PETSC_INCLUDES}
    ${HDF5_INCLUDE_DIRS}
    # ${Armadillo_INCLUDE_DIRS}
    $ENV{ARMADILLO_ROOT}/include
    ${Boost_INCLUDE_DIRS}
    ${AMGX_INCLUDE_DIRS}
)

target_link_libraries(ext_libs INTERFACE
    ${PETSC_LIBRARIES}
    ${MPI_C_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ARMADILLO_LIBRARY_PATH}
    ${Boost_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${ScaLAPACK_LIBRARIES}
    ${AMGX_LIBRARIES}    
    ${CUDA_LIBRARIES}
)

target_compile_definitions(ext_libs INTERFACE 
    ${PETSC_DEFINITIONS}
)

option(BUILD_TESTS "Build tests" ON)
option(USE_AGOS "Use AGOS for ODE solvers" OFF) 

add_subdirectory(src/util)
add_subdirectory(src/linalg)
add_subdirectory(src/mesh)
add_subdirectory(src/fem)
add_subdirectory(src/odes)
add_subdirectory(src/pdes)
add_subdirectory(src/cardiac)
add_subdirectory(src/nls)

add_subdirectory(app)

# if(BUILD_TESTS)
#    add_subdirectory(test)
# endif(BUILD_TESTS)