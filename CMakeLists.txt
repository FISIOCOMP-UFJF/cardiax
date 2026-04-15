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

# --- PETSC (via pkg-config, conda-forge compatible) ---
find_package(PkgConfig REQUIRED)

if(NOT DEFINED ENV{PKG_CONFIG_PATH})
  set(ENV{PKG_CONFIG_PATH} "$ENV{CONDA_PREFIX}/lib/pkgconfig")
else()
  set(ENV{PKG_CONFIG_PATH} "$ENV{CONDA_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

message(STATUS "PETSc found via pkg-config:")
message(STATUS "  Version  : ${PETSC_VERSION}")
message(STATUS "  Includes : ${PETSC_INCLUDE_DIRS}")
message(STATUS "  Libraries: ${PETSC_LINK_LIBRARIES}")

set(PETSC_INCLUDES    ${PETSC_INCLUDE_DIRS})
set(PETSC_LIBRARIES   ${PETSC_LINK_LIBRARIES})
set(PETSC_DEFINITIONS ${PETSC_CFLAGS_OTHER})

# --- MPI ---
find_package(MPI REQUIRED)

# --- HDF5 ---
find_package(HDF5 REQUIRED)

# --- BLAS/LAPACK (conda Netlib — must match Armadillo's build) ---
find_library(BLAS_LIBRARY
  NAMES blas
  HINTS $ENV{CONDA_PREFIX}
  PATH_SUFFIXES lib lib64
  NO_CMAKE_SYSTEM_PATH)

find_library(LAPACK_LIBRARY
  NAMES lapack
  HINTS $ENV{CONDA_PREFIX}
  PATH_SUFFIXES lib lib64
  NO_CMAKE_SYSTEM_PATH)

if(BLAS_LIBRARY AND LAPACK_LIBRARY)
  message(STATUS "BLAS   (conda): ${BLAS_LIBRARY}")
  message(STATUS "LAPACK (conda): ${LAPACK_LIBRARY}")
  set(BLAS_LIBRARIES   ${BLAS_LIBRARY})
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
else()
  message(FATAL_ERROR "conda BLAS/LAPACK not found. Run: conda install -c conda-forge blas lapack")
endif()

# --- ScaLAPACK (optional) ---
find_package(ScaLAPACK QUIET)
if(NOT ScaLAPACK_FOUND)
  set(ScaLAPACK_LIBRARIES "")
  message(STATUS "ScaLAPACK not found — linking without it.")
endif()

# --- Armadillo ---
find_library(ARMADILLO_LIBRARY
  NAMES armadillo
  HINTS $ENV{ARMADILLO_ROOT} $ENV{CONDA_PREFIX}
  PATH_SUFFIXES lib lib64)
find_path(ARMADILLO_INCLUDE_DIR armadillo
  HINTS $ENV{ARMADILLO_ROOT} $ENV{CONDA_PREFIX}
  PATH_SUFFIXES include)
if(NOT ARMADILLO_LIBRARY)
  message(FATAL_ERROR "Armadillo not found. Set ARMADILLO_ROOT or run: conda install -c conda-forge armadillo")
else()
  message(STATUS "Armadillo: ${ARMADILLO_LIBRARY}")
endif()

# # --- AmgX ---
# find_library(AMGX_LIBRARY
#   NAMES amgx
#   HINTS $ENV{AMGX_ROOT} $ENV{CONDA_PREFIX}
#   PATH_SUFFIXES lib lib64)
# find_path(AMGX_INCLUDE_DIR amgx_c.h
#   HINTS $ENV{AMGX_ROOT} $ENV{CONDA_PREFIX}
#   PATH_SUFFIXES include)
# if(NOT AMGX_LIBRARY)
#   message(FATAL_ERROR "AmgX not found. Build and install it first.")
# else()
#   message(STATUS "AmgX: ${AMGX_LIBRARY}")
# endif()

# # Eigen
# find_package(Eigen3 REQUIRED)
# if(Eigen3_FOUND)
#   message(STATUS "Eigen3 found via find_package: ${EIGEN3_INCLUDE_DIR}")
# else()
#   find_path(EIGEN3_INCLUDE_DIR Eigen/Core
#     HINTS $ENV{CONDA_PREFIX}/include/eigen3
#           $ENV{CONDA_PREFIX}/include
#     NO_CMAKE_SYSTEM_PATH)
#   if(EIGEN3_INCLUDE_DIR)
#     message(STATUS "Eigen3 found manually: ${EIGEN3_INCLUDE_DIR}")
#   else()
#     message(FATAL_ERROR "Eigen3 not found. Run: conda install -c conda-forge eigen")
#   endif()
# endif()


# # autodiff (header-only, precisa do Eigen)
# find_package(autodiff REQUIRED)
# if(autodiff_FOUND)
#   message(STATUS "autodiff found via find_package")
# else()
#   find_path(AUTODIFF_INCLUDE_DIR autodiff/forward/dual.hpp
#     HINTS $ENV{CONDA_PREFIX}/include
#     NO_CMAKE_SYSTEM_PATH)
#   if(AUTODIFF_INCLUDE_DIR)
#     message(STATUS "autodiff found manually: ${AUTODIFF_INCLUDE_DIR}")
#   else()
#     message(FATAL_ERROR "autodiff not found. Run: conda install -c conda-forge autodiff")
#   endif()
# endif()

# # --- CUDA ---
# find_library(CUDART_LIB   cudart     HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
# find_library(CUSPARSE_LIB cusparse   HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
# find_library(CUSOLVER_LIB cusolver   HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
# find_library(NVTX_LIB     nvToolsExt HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
# find_library(CUBLAS_LIB   cublas     HINTS $ENV{CUDAToolkit_ROOT} PATH_SUFFIXES lib lib64)
# set(CUDA_LIBRARIES ${CUDART_LIB} ${CUSPARSE_LIB} ${CUSOLVER_LIB} ${NVTX_LIB} ${CUBLAS_LIB})

# --- Interface library ---
add_library(ext_libs INTERFACE)

target_include_directories(ext_libs INTERFACE
    ${PROJECT_SOURCE_DIR}/src
    ${MPI_C_INCLUDE_DIRS}
    ${PETSC_INCLUDES}
    ${HDF5_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIR}
    $ENV{ARMADILLO_ROOT}/include
    # ${AMGX_INCLUDE_DIR}
    # ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(ext_libs INTERFACE
    ${PETSC_LIBRARIES}
    ${MPI_C_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ARMADILLO_LIBRARY}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${ScaLAPACK_LIBRARIES}
    # ${AMGX_LIBRARY}
    # ${CUDA_LIBRARIES}
)

target_compile_definitions(ext_libs INTERFACE
    ${PETSC_DEFINITIONS}
)

# --- Subdirectories ---
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
# endif()
