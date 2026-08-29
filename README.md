# C A R D I A X

This document provides installation instructions and guidance for configuring the environment required to build and run Cardiax.

## Contents

1. [Prerequisites](#1-prerequisites)
2. [Installation](#2-installation)
   - [2.1 Automatic Installation With Conda](#21-automatic-installation-with-conda)
   - [2.2 Manual Installation With Conda](#22-manual-installation-with-conda)
   - [2.3 Manual Compilation of Libraries](#23-manual-compilation-of-libraries)
3. [Documentation](#3-documentation)
4. [References](#4-references)

---

## 1. [Prerequisites](#1-prerequisites)

Before starting, make sure the following dependencies are installed on your system:

* **Miniconda** or **Anaconda**
* **Nvidia Drivers** compatible with CUDA 12
* **Git**

Dependencies:
*  **PETSc**
*  **Armadillo**
*  **HDF5**
*  **AMGx (optional, for GPU solvers)**

--- 

## 2. [Installation](#2-installation)

- [2.1 Automatic Installation With Conda](#21-automatic-installation-with-conda)
- [2.2 Manual Installation With Conda](#22-manual-installation-with-conda)
- [2.3 Manual Compilation of Libraries](#23-manual-compilation-of-libraries)

---

### [2.1 Automatic Installation With Conda](#21-automatic-installation-with-conda)

Cardiax uses a set of automated scripts to handle complex dependencies (HDF5, PETSc, AMGX and Armadillo) and the compilation process within an isolated Conda environment.

**0. Deactivate Conda Environments**
Before running the installation scripts, it is crucial to exit any active Conda environment to avoid conflicts. Run the following command (multiple times if necessary) until you are entirely out of any Conda environment:

```
conda deactivate
```

1. Grant execution permissions to the .sh scripts in the project root:

```
chmod +x ./install_deps.sh ./build.sh
```

2. Run the dependency installation script. This step usually needs to be performed only once. It will download, compile and configure all dependencies inside a temporary Conda environment.

By default, running the script without arguments installs all dependencies:

```
./install_deps.sh 
```
*Note: This step may take some time, as several libraries are compiled from source — so grab a coffee.*

Custom Installation: If you want to install or skip specific dependencies individually, you can use the force and ignore flags (e.g., --force-hdf5, --ignore-amgx). To see all available configuration flags, run:  

```
./install_deps.sh --help
```

3. Once the dependencies are installed, compile the simulator using the build script: 
```
./build.sh
```
*If you make changes to the source code, you only need to rerun this command to recompile the Cardiax binary.*

You can also specify a target binary name to compile only that executable. For example:

```
./build.sh electromech
```
*To list all available executables, use the `--help` flag.*

If you have AMGX installed but wish to bypass it and compile using only the PETSc CPU solver, use the --force-petsc-solver flag:
```
./build.sh --force-petsc-solver
```
*(Note: If AMGX is not installed or cannot be found on your system, the build process will automatically fall back to the PETSc solver)*


4. After a successful build, all executables will be available in the build/app directory.
For more detailed usage examples, refer to the README located in the examples/ folder.
If you just want to get things started, you can run:

```
./build/app/electromech -f examples/pvloop.xml -amgx configs/CG_DILU.json
```

#### Dependency Management & Troubleshooting
1. The install_deps.sh script performs checks to determine whether each dependency has already been compiled, avoiding unnecessary rebuilds. However, if you need to update or force the recompilation of specific dependencies, several force flags are available.

```
./install_deps.sh --help
```

2. By default, install_deps.sh and the corresponding conda recipes install AMGX targeting a specific set of GPU architectures (70, 75, 80, 86, 89).
If your GPU architecture is not included in the conda recipe, the compilation will succeed, but you may encounter runtime issues when executing Cardiax.
To fix this, locate the CUDA_ARCH variable in recipe/amgx_custom/build.sh and add the required architecture.

---

### [2.2 Manual Installation With Conda](#22-manual-installation-with-conda)

This installation method allows you to manually create the Conda environment and install the required libraries, while leaving the compilation of Cardiax to CMake.

First, create and activate a dedicated Conda environment:

```bash
conda create -c conda-forge -n cardiax-env
conda activate cardiax-env
```

Install the required dependencies from `conda-forge`:

```bash
conda install conda-forge::petsc
conda install conda-forge::armadillo
conda install conda-forge::hdf5
```

Create the build directory and configure the project using CMake:

```bash
mkdir build
cd build
cmake ..
```

During the CMake configuration step, verify that the libraries are detected from the active Conda environment. A successful configuration should produce output similar to:

```text
-- Checking for module 'PETSc'
--   Found PETSc, version 3.25.4
-- PETSc found via pkg-config:
--   Version  : 3.25.4
--   Includes : /home/rocha/miniconda3/envs/cardiax-env/include
--   Libraries: /home/rocha/miniconda3/envs/cardiax-env/lib/libpetsc.so
-- Found MPI_C: /home/rocha/miniconda3/envs/cardiax-env/lib/libmpi.so (found version "5.0")
-- Found MPI_CXX: /home/rocha/miniconda3/envs/cardiax-env/lib/libmpicxx.so (found version "5.0")
-- Found MPI: TRUE (found version "5.0")
-- Found HDF5: hdf5-shared (found version "2.2.0")
-- BLAS   (conda): /home/rocha/miniconda3/envs/cardiax-env/lib/libblas.so
-- LAPACK (conda): /home/rocha/miniconda3/envs/cardiax-env/lib/liblapack.so
-- Armadillo: /home/rocha/miniconda3/envs/cardiax-env/lib/libarmadillo.so
-- AmgX not found. Falling back to PETSC_SOLVER.
```

> **Note:** The paths shown above are examples and will depend on where your Conda or Miniconda installation is located.

If all dependencies are detected successfully, compile Cardiax with:

```bash
make
```

or 

```bash
make -j8
```

AMGX is optional. If it is not installed or cannot be found, Cardiax automatically falls back to the PETSc CPU solver.

---

### [2.3 Manual Compilation of Libraries](#23-manual-compilation-of-libraries)

If you are deploying Cardiax on an HPC cluster or prefer not to use Conda, you can build all dependencies manually from source. The source files for the specific versions of Armadillo, HDF5, and PETSc are located in the `dependencies/`folder. 

**0. Setup Installation Directory**
First, define an installation directory (prefix) where all libraries will be installed, and set up your compiler variables.  
```
export PREFIX="/path/to/your/custom/install/dir"
export CC=gcc
export CXX=g++
export FC=gfortran
```

***1. Install Armadillo***

```
cd dependencies
tar -xzf armadillo-*.tar.gz
cd armadillo-*
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX
make 
make install
cd ../..
```

***2. Install HDF5***
```
cd dependencies
tar -xzf hdf5-*.tar.gz
cd hdf5-*
./configure --prefix=$PREFIX
make 
make install
cd ../..
```

***3. Install PETSc***
Extract PETSc and run the configure script. It will automatically download and build internal dependencies like MPICH, MUMPS and HYPRE: 

```
cd dependencies
tar -xzf petsc-*.tar.gz
cd petsc-*

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
  LIBS="-lgfortran -lquadmath -lgomp"

make PETSC_DIR=$PWD PETSC_ARCH=arch-linux2-c-debug all
make PETSC_DIR=$PWD PETSC_ARCH=arch-linux2-c-debug install
cd ../..
```


***4. Install AMGX (Optional but Recommended)***
If you want to use GPU-accelerated linear solvers, you must download and compile Nvidia's AMGX from their official repository:

```
cd dependencies
git clone [https://github.com/NVIDIA/AMGX.git](https://github.com/NVIDIA/AMGX.git)
cd AMGX
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX
make 
make install
cd ../..
```

***5. Building Cardiax***
Since you are not using Conda, do not use the build.sh script, as it assumes a Conda environment. Instead, create a build directory in the Cardiax root and run CMake manually, pointing it to your newly installed libraries:

```
mkdir build && cd build

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DPETSC_DIR=$PREFIX \
  -DAMGX_DIR=$PREFIX \
  -DARMADILLO_DIR=$PREFIX \
  -DHDF5_ROOT=$PREFIX \
  -DCUDA_DIR=/cuda12/directory

make 
```

---

## 3. [Documentation](#3-documentation)

Code documentation and examples will be provided soon.

---

## 4. [References](#4-references)

Reference list:

- Campos, Joventino Oliveira, et al. "Preconditioned augmented Lagrangian formulation for nearly incompressible cardiac mechanics." *International Journal for Numerical Methods in Biomedical Engineering* 34.4 (2018): e2948. [Link](https://onlinelibrary.wiley.com/doi/abs/10.1002/cnm.2948)
- Uncertainty quantification and sensitivity analysis of left ventricular function during the full cardiac cycle. [Link](https://royalsocietypublishing.org/doi/abs/10.1098/rsta.2019.0381)
- Effects of left ventricle wall thickness uncertainties on cardiac mechanics.* [Link](https://link.springer.com/article/10.1007/s10237-019-01153-1)
- Polynomial chaos expansion surrogate modeling of passive cardiac mechanics using the Holzapfel–Ogden constitutive model. [Link](https://www.sciencedirect.com/science/article/abs/pii/S1877750323000996)

