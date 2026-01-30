#!/bin/bash
set -e

echo ">>> COMPILANDO AMGX (CUSTOM) <<<"

cd $SRC_DIR
rm -rf build
mkdir -p build && cd build

# Configuração Universal (Fat Binary)
# Removemos -DCMAKE_CUDA_ARCHITECTURES para ele compilar para todas as placas suportadas.
# Desligamos NVTX e MPI interno para evitar conflitos.

CUDA_ARCH="${CARDIAX_CUDA_ARCH:-70;75;80;86;89}"

cmake .. \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_CUDA_FLAGS="-ccbin $CXX -allow-unsupported-compiler" \
    -DAMGX_BUILD_TESTS=OFF \
    -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE \
    -DNVTX=OFF \
    -DENABLE_NVTX=OFF \
    -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCH"

make -j$CPU_COUNT
make install
