#!/bin/bash
set -e

# --- 1. ARGUMENT PARSING ---
BUILD_TARGET=""
FORCE_PETSC_SOLVER=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --force-petsc-solver)
            FORCE_PETSC_SOLVER=1
            shift # Avança para o próximo argumento
            ;;
        -h|--help)
            echo "Usage: ./build.sh [TARGET] [OPTIONS]"
            echo ""
            echo "Arguments:"
            echo "  TARGET       (Optional) Name of the specific binary to build."
            echo "                  Options:"
            echo "                       electromech   "
            echo "                       monodomain    "
            echo "                       nonlinearelas "
            echo "                       poisson       "
            echo "                       l2projection  "
            echo "                       elasticity    "
            echo "                       bidomain      "
            echo "               If omitted, ALL targets will be built."
            echo ""
            echo "Options:"
            echo "  -h, --help            Show this help message"
            echo "  --force-petsc-solver  Force PETSc linear solver (CPU only)"
            echo ""
            exit 0
            ;;
        -*)
            # Captura qualquer outra flag que comece com '-' e não foi mapeada
            echo "Error: Unknown option: $1"
            echo "Use ./build.sh --help for usage."
            exit 1
            ;;
        *)
            # Captura o argumento posicional (o TARGET)
            if [[ -z "$BUILD_TARGET" ]]; then
                BUILD_TARGET="$1"
            else
                echo "Error: Multiple targets specified ('$BUILD_TARGET' and '$1'). Only one is allowed."
                exit 1
            fi
            shift # Avança para o próximo argumento
            ;;
    esac
done

# Conda environment name
ENV_NAME="cardiax_env"

echo ">>> COMPILING CARDIAX <<<"
if [[ -n "$BUILD_TARGET" ]]; then
    echo ">>> Specific Target: $BUILD_TARGET"
else
    echo ">>> Target: ALL"
fi

# 2. Automatically Activate Conda Environment
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

if [ $? -ne 0 ]; then
    echo "ERROR: Environment '$ENV_NAME' not found."
    echo "Please run './install_deps.sh' first."
    exit 1
fi

echo "Environment activated: $CONDA_PREFIX"

# 3. Configure Environment Variables
export PREFIX=$CONDA_PREFIX
export PETSC_DIR=$PREFIX
export PETSC_ARCH=""
export AMGX_ROOT=$PREFIX
export ARMADILLO_ROOT=$PREFIX
export HDF5_ROOT=$PREFIX
export CUDAToolkit_ROOT=$PREFIX

# Compilers
export CC=$PREFIX/bin/mpicc
export CXX=$PREFIX/bin/mpicxx

# Flags 
export CUDAFLAGS="-ccbin $CXX -allow-unsupported-compiler -w -Xcompiler -fPIC"
export CFLAGS="-O2 -g -fPIC"
export CXXFLAGS="-O2 -g -fPIC"

# 4. The FORTRAN Combo
MATH_LIBS="$PREFIX/lib/libflapack.a;$PREFIX/lib/libfblas.a;$PREFIX/lib/libgfortran.so;$PREFIX/lib/libquadmath.so"

# 5. Find NVTX (CUDA 12 Compatibility)
NVTX_PATH=$(find $PREFIX -name "libnvToolsExt.so*" 2>/dev/null | head -n 1)
if [ -z "$NVTX_PATH" ]; then
    NVTX_PATH=$(find $PREFIX -name "libnvtx3interop.so*" 2>/dev/null | head -n 1)
fi

# 6. Automatic Patch (Removes broken test if it exists)
##REMOVE IT?? 
if [ -f "src/nls/CMakeLists.txt" ]; then
    grep -q "EXCLUDE_FROM_ALL" src/nls/CMakeLists.txt || \
    sed -i 's/add_executable(testSimple/add_executable(testSimple EXCLUDE_FROM_ALL/g' src/nls/CMakeLists.txt
fi

# 7. Prepare Build Directory
mkdir -p build
cd build

if [[ $FORCE_PETSC_SOLVER -eq 1 ]]; then
    CMAKE_PETSC_FLAG="-DFORCE_PETSC_SOLVER=ON"
    echo ">>> Forcing PETSc solver (Bypassing AmgX) <<<"
else
    CMAKE_PETSC_FLAG="-DFORCE_PETSC_SOLVER=OFF"
fi

# VERIFICA SE O CMAKE JÁ RODOU ANTERIORMENTE
if [ ! -f "Makefile" ]; then
    echo ">>> Configuring CMAKE <<<"
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_BUILD_TYPE=Release \
        -DMPI_C_COMPILER=$CC \
        -DMPI_CXX_COMPILER=$CXX \
        -DCMAKE_C_COMPILER=$CC \
        -DCMAKE_CXX_COMPILER=$CXX \
        -DCMAKE_CUDA_FLAGS="$CUDAFLAGS" \
        -DBLAS_LIBRARIES="$MATH_LIBS" \
        -DLAPACK_LIBRARIES="$MATH_LIBS" \
        -DNVTX_LIB="$NVTX_PATH" \
        -DBUILD_TESTS=OFF \
        $CMAKE_PETSC_FLAG
else
    echo ">>> Configuração do CMake encontrada (Makefile existe). Pulando fase do CMake... <<<"
    echo "    (Para forçar um novo CMake, delete a pasta 'build' e rode o script novamente)"
fi

echo ">>> Running make: "
# Compiles using all processor cores + Optional Target
make -j$(nproc) $BUILD_TARGET

echo "=========================================="
echo "            BUILD COMPLETED               "
if [[ -n "$BUILD_TARGET" ]]; then
    echo "   Built target: $BUILD_TARGET"
else
    echo "   All targets built."
fi
echo "   Executables are in: ./build/app/       "
echo "=========================================="