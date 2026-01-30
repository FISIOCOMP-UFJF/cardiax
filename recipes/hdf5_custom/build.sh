#!/bin/bash
set -e

echo ">>> COMPILANDO HDF5 1.8.13 (CUSTOM) <<<"

# Flags cruciais para GCC 12+
export CFLAGS="-O2 -std=gnu99 -fcommon -w -Wno-incompatible-pointer-types -Wno-int-conversion -Wno-implicit-function-declaration -Wno-implicit-int"

# Configuração
./configure \
    --prefix=$PREFIX \
    --enable-shared \
    --enable-static=no \
    --enable-production \
    --disable-debug \
    --with-zlib=no \
    CFLAGS="$CFLAGS"

# Compilação e Instalação
make -j$CPU_COUNT
make install

# Limpeza opcional para economizar espaço no pacote final
rm -rf $PREFIX/share/hdf5_examples
