#!/bin/bash
set -e

echo ">>> COMPILANDO ARMADILLO (CUSTOM) <<<"

# O Conda já descompacta o código na pasta $SRC_DIR
cd $SRC_DIR

# Configuração do CMake
cmake . \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5

# Compilação e Instalação
make -j$CPU_COUNT
make install
