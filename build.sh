#!/bin/bash
set -e

# --- 1. ARGUMENT PARSING ---
BUILD_TARGET=""

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: ./build.sh [TARGET]"
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
    exit 0
elif [[ -n "$1" ]]; then
    BUILD_TARGET="$1"
fi

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

echo ">>> Configuring CMAKE <<<"
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE:-Release} \
    -DMPI_C_COMPILER=$CC \
    -DMPI_CXX_COMPILER=$CXX \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_CUDA_FLAGS="$CUDAFLAGS" \
    -DBLAS_LIBRARIES="$MATH_LIBS" \
    -DLAPACK_LIBRARIES="$MATH_LIBS" \
    -DNVTX_LIB="$NVTX_PATH" \
    -DBUILD_TESTS=OFF

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