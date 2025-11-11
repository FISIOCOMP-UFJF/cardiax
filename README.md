# C A R D I A X

This document provides installation instructions and guidance for configuring the environment required to build and run CardiaX.
---

## 1. Requirements

CardiaX depends on the following external libraries:

- HDF5 **1.8.13**
- Armadillo (compatible with C++11 or later)
- AMGx (NVIDIA CUDA-based algebraic multigrid)
- PETSc **3.19.\*** (this guide uses **3.19.6** as an example)

### General Recommendations

1. Create a directory to store all external dependencies, for example:
```
mkdir -p /home/<username>/source
```
2. Always export environment variables in your `~/.bashrc`.
3. For Armadillo, make sure to export the path **containing both `include/` and `lib/`** directories. Depending on your system, this folder might be named `armadillo`, `usr`, or `usr/local`.

---

## 2. Installing Dependencies

### 2.1 HDF5 (1.8.13)

```bash
tar xvf hdf5-1.8.13.tar
cd hdf5-1.8.13
./configure --prefix=/home/<username>/source/hdf5
make -j4
make install
```

Add to your ~/.bashrc:

```
export HDF5_ROOT=/home/<username>/source/hdf5
```

### 2.2 AMGX
```
git clone --recursive https://github.com/NVIDIA/AMGX
mv AMGX amgx
cd amgx
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/<username>/source/amgx \
      -DCMAKE_C_COMPILER=<path_to_gcc> \
      -DCMAKE_CXX_COMPILER=<path_to_g++> \
      ..
make -j16 all install
```
Add to your ~/.bashrc
```
export AMGX_ROOT=/home/<username>/source/amgx
```


### 2.3 Armadillo

```
tar -zxf armadillo-8.300.1.tar.gz
cd armadillo-8.300.1
mkdir build && cd build
cmake ..
make -j4
make install DESTDIR=/home/<username>/source/armadillo
```

Add to your ~/.bashrc
(ensure this path contains the include/ and lib/ folders):

```
export ARMADILLO_ROOT=/home/<username>/source/armadillo/usr
```

### 2.4 PETSc (3.19.*)

```
tar xvf petsc-3.19.6.tar.gz
cd petsc-3.19.6/

python3 ./configure \
 --COPTFLAGS="-O2 -g" --CXXOPTFLAGS="-O2 -g" --FOPTFLAGS="-O2 -g -std=legacy" \
 --with-debugging=yes \
 --download-mpich \
 --download-suitesparse \
 --with-hdf5-dir=$HDF5_ROOT \
 --download-fblaslapack \
 --download-hypre \
 --download-mumps \
 --download-scalapack \
 --download-superlu

make PETSC_DIR=/home/<username>/source/petsc-3.19.6 PETSC_ARCH=arch-linux2-c-debug all
```

Add to your ~/.bashrc:

```
export PETSC_DIR=/home/<username>/source/petsc-3.19.6
export PETSC_ARCH=arch-linux2-c-debug
```

## 3. Building CardiaX

```
git clone https://github.com/FISIOCOMP-UFJF/cardiax.git
cd cardiax
mkdir build && cd build
cmake ..
make -j4
```


## 4.PETSc tips

### 3D elasticity preconditioner command line options

#### BoomerAMG
```
-ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_iter 1 -ksp_view -ksp_monitor
```

#### GAMG
```
-ksp_type cg -pc_type gamg -pc_gamg_type agg -log_summary
-ksp_monitor -ksp_view -options_left -mg_levels_ksp_max_it 1
```
ou
```
-ksp_type cg -pc_type gamg -pc_gamg_type agg -log_summary
-ksp_monitor -ksp_view -options_left -mg_levels_ksp_type richardson -mg_levels_pc_type sor
```




