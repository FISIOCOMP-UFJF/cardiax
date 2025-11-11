# CARDIAX
# author: Bernardo Martins Rocha
# mail: bernardomartinsrocha@gmail.com

cmake_minimum_required(VERSION 3.12) # 3.10 funciona, mas 3.12+ é melhor para CMake moderno
project(Cardiax C CXX)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 OLD)
endif()

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()

# --- Configurações Iniciais ---
# Remove flags globais. Deixe o CMake gerenciar isso com Build Types.
# O usuário pode compilar com: cmake -DCMAKE_BUILD_TYPE=Release ..
# Ou -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_CXX_STANDARD 17) # Exemplo: defina um padrão C++
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configura cmake modules dir
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

# --- Configuração de Pacotes Externos ---
# Esta seção está boa! find_package é a maneira correta.
if (NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "A variável de ambiente PETSC_DIR deve estar configurada para o PETSc.")
endif()
if (NOT DEFINED ENV{PETSC_ARCH})
    message(WARNING "A variável de ambiente PETSC_ARCH não está configurada. Usando 'arch-linux2-c-debug' como padrão.")
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

# Adiciona includes das dependências
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

# Adiciona bibliotecas das dependências
target_link_libraries(ext_libs INTERFACE
    ${PETSC_LIBRARIES}
    ${MPI_C_LIBRARIES}
    ${HDF5_LIBRARIES}
    $ENV{ARMADILLO_ROOT}/lib64/libarmadillo.so
    ${Boost_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${ScaLAPACK_LIBRARIES}
    ${AMGX_LIBRARIES}    
    ${CUDA_LIBRARIES}
)

# Adiciona definições das dependências
target_compile_definitions(ext_libs INTERFACE 
    ${PETSC_DEFINITIONS}
)


# --- [REMOVIDO] ---
# O `include_directories(...)` global foi removido.
# O `set(EXTERNAL_LIBRARIES ...)` global foi removido.
# O `add_definitions(${PETSC_DEFINITIONS})` global foi removido.
# Tudo isso agora está dentro do alvo "ext_libs".


# --- Subdiretórios e Alvos ---
option(BUILD_TESTS "Build tests" ON)
option(USE_AGOS "Use AGOS for ODE solvers" OFF) # Você pode usar isso para CFLAGS condicionais

# Adiciona os subdiretórios que definem seus próprios alvos (bibliotecas)
add_subdirectory(src/util)
add_subdirectory(src/linalg)
add_subdirectory(src/mesh)
add_subdirectory(src/fem)
add_subdirectory(src/odes)
add_subdirectory(src/pdes)
add_subdirectory(src/cardiac)
add_subdirectory(src/nls)

# O aplicativo principal, que depende das bibliotecas acima
add_subdirectory(app)

# if(BUILD_TESTS)
#    add_subdirectory(test)
# endif(BUILD_TESTS)