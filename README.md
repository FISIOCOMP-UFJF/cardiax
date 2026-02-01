# C A R D I A X

This document provides installation instructions and guidance for configuring the environment required to build and run Cardiax.

---

## 1. Prerequisites

Before starting, make sure the following dependencies are installed on your system:

* **Miniconda** or **Anaconda**
* **Nvidia Drivers** compatible with CUDA 12
* **Git**

--- 

## 2. Instalation 
Cardiax uses a set of automated scripts to handle complex dependencies (HDF5, PETSc, AMGX and Armadillo) and the compilation process within an isolated Conda environment.

1. Grant execution permissions to the .sh scripts in the project root:

```
chmod +x ./install_deps.sh ./build.sh
```

2. Run the dependency installation script. This step usually needs to be performed only once. It will download, compile and configure all dependencies inside a temporary Conda environment.

```
./install_deps.sh 
```
*Note: This step may take some time, as several libraries are compiled from source — so grab a coffee.*

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


4. After a successful build, all executables will be available in the build/app directory.
For more detailed usage examples, refer to the README located in the examples/ folder.
If you just want to get things started, you can run:

```
./build/app/electromech -f examples/pvloop.xml -s ul -amgx configs/CG_DILU.json
```

## 3. Dependency Management & Troubleshooting
1. The install_deps.sh script performs checks to determine whether each dependency has already been compiled, avoiding unnecessary rebuilds. However, if you need to update or force the recompilation of specific dependencies, several force flags are available.

```
./install_deps.sh --help
```

2. By default, install_deps.sh and the corresponding conda recipes install AMGX targeting a specific set of GPU architectures (70, 75, 80, 86, 89).
If your GPU architecture is not included in the conda recipe, the compilation will succeed, but you may encounter runtime issues when executing Cardiax.
To fix this, locate the CUDA_ARCH variable in recipe/amgx_custom/build.sh and add the required architecture.