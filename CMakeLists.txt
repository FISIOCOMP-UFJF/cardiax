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

# --- CONFIGURAÇÃO PETSC ---
if (NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "The environment variable PETSC_DIR must be set to the PETSc directory.")
endif()
if (NOT DEFINED ENV{PETSC_ARCH})
    message(WARNING "The environment variable PETSC_ARCH is not set. Using 'arch-linux2-c-debug' as the default.")
    set(ENV{PETSC_ARCH} "arch-linux2-c-debug")
endif()

# Caminho onde o PETSc instala dependências (MPI, BLAS, LAPACK, etc)
set(PETSC_DEPS_ROOT "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")

find_package(PETSc REQUIRED)

# --- BUSCA DE DEPENDÊNCIAS COM FALLBACK PARA PETSC ---

# 1. MPI
find_package(MPI QUIET)
if(NOT MPI_FOUND)
    message(STATUS "System MPI not found. Searching in PETSc dir: ${PETSC_DEPS_ROOT}")
    set(ENV{MPI_HOME} "${PETSC_DEPS_ROOT}")
    list(APPEND CMAKE_PREFIX_PATH "${PETSC_DEPS_ROOT}")
    find_package(MPI REQUIRED)
endif()

find_package(HDF5 REQUIRED)

# 2. BLAS
find_package(BLAS QUIET)
if(NOT BLAS_FOUND)
    message(STATUS "System BLAS not found. Searching in PETSc dir...")
    list(APPEND CMAKE_PREFIX_PATH "${PETSC_DEPS_ROOT}")
    find_package(BLAS REQUIRED)
endif()

# 3. LAPACK
find_package(LAPACK QUIET)
if(NOT LAPACK_FOUND)
    message(STATUS "System LAPACK not found. Searching in PETSc dir...")
    list(APPEND CMAKE_PREFIX_PATH "${PETSC_DEPS_ROOT}")
    find_package(LAPACK REQUIRED)
endif()

# 4. ScaLAPACK
find_package(ScaLAPACK QUIET)
if(NOT ScaLAPACK_FOUND)
    message(STATUS "System ScaLAPACK not found. Searching in PETSc dir...")
    list(APPEND CMAKE_PREFIX_PATH "${PETSC_DEPS_ROOT}")
    find_package(ScaLAPACK REQUIRED)
endif()

# --- OUTRAS DEPENDÊNCIAS ---
# find_package(Armadillo REQUIRED)

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
    ${AMGX_INCLUDE_DIRS}
)

target_link_libraries(ext_libs INTERFACE
    ${PETSC_LIBRARIES}
    ${MPI_C_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ARMADILLO_LIBRARY_PATH}
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
#   add_subdirectory(test)
# endif(BUILD_TESTS)