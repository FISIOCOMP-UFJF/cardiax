#!/bin/bash
set -e

echo ">>> COMPILANDO PETSC 3.19.6 (CUSTOM) <<<"

# Garante que o compilador ache o HDF5
export CPATH=$PREFIX/include:$CPATH
export LIBRARY_PATH=$PREFIX/lib:$LIBRARY_PATH

# Configuração (Seguindo seu README)
# Python 3.10 (do ambiente) roda o configure sem erro xdrlib
python3 ./configure \
  --prefix=$PREFIX \
  --with-petsc-arch=arch-linux2-c-debug \
  --COPTFLAGS="-O2 -g" \
  --CXXOPTFLAGS="-O2 -g" \
  --FOPTFLAGS="-O2 -g -std=legacy" \
  --with-debugging=yes \
  --download-mpich \
  --download-fblaslapack \
  --download-hypre \
  --download-suitesparse \
  --download-mumps \
  --download-cmake \
  --download-scalapack \
  --download-superlu \
  --with-hdf5-dir=$PREFIX \
  --with-fc=$FC --with-cc=$CC --with-cxx=$CXX \
  --with-make-np=$CPU_COUNT \
  LIBS="-lgfortran -lquadmath -lgomp"

# Compilação
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux2-c-debug all

# Instalação
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux2-c-debug install
