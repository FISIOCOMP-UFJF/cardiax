#!/bin/bash
set -e

echo ">>> COMPILANDO PETSC 3.19.6 (CUSTOM) <<<"

# Garante que o compilador ache o HDF5
export CPATH=$PREFIX/include:$CPATH
export LIBRARY_PATH=$PREFIX/lib:$LIBRARY_PATH

# Variante CUDA + Kokkos (padrao do install_deps.sh; --petsc-no-kokkos usa a variante cpu)
# Habilita VECKOKKOS / MATAIJKOKKOS (-vec_type kokkos -mat_type aijkokkos).
# O Kokkos e o Kokkos Kernels ficam instalados em $PREFIX (lib/cmake/Kokkos), entao
# o proprio Cardiax pode usar o MESMO Kokkos via find_package(Kokkos).
KOKKOS_OPTS=()
MAKE_NP=$CPU_COUNT
# Variante cpu: MPICH e hypre baixados pelo PETSc, como sempre foi.
MPI_OPTS=(--download-mpich --with-fc=$FC --with-cc=$CC --with-cxx=$CXX)
HYPRE_OPTS=(--download-hypre)
if [ "${CARDIAX_PETSC_KOKKOS:-0}" = "1" ]; then
  ARCH="${CARDIAX_PETSC_CUDA_ARCH:?CARDIAX_PETSC_CUDA_ARCH nao definido}"
  echo ">>> Variante CUDA + Kokkos, sm_${ARCH} <<<"
  # Layout do CUDA do conda-forge: headers em targets/x86_64-linux/include, libs em lib/,
  # stub do libcuda em lib/stubs. Libs passadas explicitamente porque o CUDA 12.9 nao tem
  # mais libnvToolsExt (NVTX3 e so header, em include/nvtx3), que o configure 3.19 procura.
  CUDA_INC="$PREFIX/targets/x86_64-linux/include"

  # O configure 3.19 acha a raiz do CUDA procurando <dir do nvcc>/../include/cuda.h, mas no
  # conda-build o nvcc fica no env de build e os headers no env host. Nao da para usar
  # --with-cuda-dir junto com --with-cuda-lib, entao montamos uma raiz "falsa" com os dois.
  # Ela tambem e passada adiante (CUDA_HOME etc.) para o CMake do Kokkos/Kokkos Kernels.
  CUDA_ROOT="$SRC_DIR/_cuda_root"
  rm -rf "$CUDA_ROOT"
  mkdir -p "$CUDA_ROOT/bin"
  # wrapper (nao symlink): o nvcc acha cicc/nvvm a partir do caminho pelo qual foi chamado.
  # -Wno-deprecated-gpu-targets: o nvcc 12.9 avisa no stderr sobre a arch padrao (< sm_75) em
  # toda chamada, e o configure 3.19 trata qualquer stderr do preprocessador como falha.
  printf '#!/bin/bash\nexec "%s" -Wno-deprecated-gpu-targets "$@"\n' "$(command -v nvcc)" > "$CUDA_ROOT/bin/nvcc"
  chmod +x "$CUDA_ROOT/bin/nvcc"
  ln -s "$CUDA_INC" "$CUDA_ROOT/include"
  ln -s "$PREFIX/lib" "$CUDA_ROOT/lib64"
  ln -s "$PREFIX/lib" "$CUDA_ROOT/lib"
  export CUDA_HOME="$CUDA_ROOT" CUDA_PATH="$CUDA_ROOT" CUDAToolkit_ROOT="$CUDA_ROOT"

  # Jobs paralelos limitados pela RAM: cada nvcc do Kokkos/Kokkos Kernels pode passar de 2-3 GB.
  # Com 12 jobs numa maquina de 15 GB o build travou o PC (swap). Usa ~1 job por 3 GB de RAM.
  # Pode ser sobrescrito com CARDIAX_PETSC_MAKE_NP.
  MEM_GB=$(awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo)
  MAKE_NP=$(( MEM_GB / 3 )); [ "$MAKE_NP" -lt 1 ] && MAKE_NP=1
  [ "$MAKE_NP" -gt "$CPU_COUNT" ] && MAKE_NP=$CPU_COUNT
  MAKE_NP=${CARDIAX_PETSC_MAKE_NP:-$MAKE_NP}
  echo ">>> make-np = $MAKE_NP (RAM ${MEM_GB} GB) <<<"

  # MPICH: com CUDA, o --download-mpich do 3.19 recompila o MPICH "CUDA-aware" sobre UCX
  # (ch4:ucx). Usamos o mpich do conda-forge (que o cardiax_env ja instala de qualquer jeito).
  MPI_OPTS=(--with-cc="$PREFIX/bin/mpicc" --with-cxx="$PREFIX/bin/mpicxx" --with-fc="$PREFIX/bin/mpifort")
  # hypre: o hypre.py do 3.19 aborta com CUDA >= 12 ("cannot be used with CUDA 12 or higher").
  # O Cardiax nao usa hypre no codigo (so via -pc_type hypre na linha de comando). AMG no
  # PETSc fica com o GAMG, que funciona com aijkokkos.
  HYPRE_OPTS=()

  # CUDA >= 12.6 (CCCL 2.x): thrust::tuple virou cuda::std::tuple, que nao tem .get<N>().
  # O PETSc 3.19 usa t.get<N>() em alguns .cu/.hpp; trocamos por thrust::get<N>(t), que
  # funciona com Thrust antigo e novo (versoes novas do PETSc ja fazem isso).
  perl -pi -e 's/\b([A-Za-z_]\w*(?:\.get_iterator_tuple\(\))?)\.get<(\d)>\(\)/thrust::get<$2>($1)/g' \
    src/mat/impls/aij/seq/seqcusparse/aijcusparse.cu \
    src/mat/impls/aij/mpi/mpicusparse/mpiaijcusparse.cu \
    src/vec/vec/impls/seq/cupm/vecseqcupm.hpp \
    src/sys/classes/random/impls/curand/curand2.cu

  KOKKOS_OPTS=(
    --with-cuda=1
    --with-cudac="$CUDA_ROOT/bin/nvcc"
    --with-cuda-arch="$ARCH"
    --with-cuda-include="[$CUDA_INC,$CUDA_INC/nvtx3]"
    --with-cuda-lib="-L$PREFIX/lib -lcudart -lcufft -lcublas -lcusparse -lcusolver -lcurand -L$PREFIX/lib/stubs -lcuda"
    --download-kokkos
    --download-kokkos-kernels
    # o libkokkos*.so depende do libmpicxx do mpich (dependencia indireta); sem rpath-link
    # o ld so avisa, mas o configure 3.19 trata o aviso como falha ("could not be used")
    LDFLAGS="-Wl,-rpath-link,$PREFIX/lib"
  )
fi

# Configuração (Seguindo seu README)
# Python 3.10 (do ambiente) roda o configure sem erro xdrlib
python3 ./configure \
  --prefix=$PREFIX \
  --with-petsc-arch=arch-linux2-c-opt \
  --COPTFLAGS="-O3" \
  --CXXOPTFLAGS="-O3" \
  --FOPTFLAGS="-O3 -std=legacy" \
  --with-debugging=0 \
  "${MPI_OPTS[@]}" \
  --download-fblaslapack \
  "${HYPRE_OPTS[@]}" \
  --download-suitesparse \
  --download-mumps \
  --download-cmake \
  --download-scalapack \
  --download-superlu \
  --with-hdf5-dir=$PREFIX \
  --with-make-np=$MAKE_NP \
  "${KOKKOS_OPTS[@]}" \
  LIBS="-lgfortran -lquadmath -lgomp"

# Compilação
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux2-c-opt all

# Instalação
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux2-c-opt install
