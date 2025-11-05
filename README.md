# C A R D I A X

## 1. Installation

     * Armadillo, BLAS and LAPACK
     * MPI (for use with PETSc)
     * PETSc

### 1.2. PETSc

You may need MPI installed to compile PETSc.

	tar xvf petsc-3.4.3.tar.gz
	cd petsc-3.4.3/

To compile PETSc WITHOUT MPI
```
python2 './configure' '--with-cc=gcc' '--with-cxx=g++'
                              '--with-fc=gfortran'
                              '--download-fblaslapack'
                              '--download-suitesparse'
                              '--with-mpi=0'
                              '--with-debugging=1'
```                              

To compile PETSc WITH MPI you can either use the following
command line and let the configure find you MPI automatically
```
python2 './configure' '--with-cc=gcc' '--with-cxx=g++'
                              '--with-fc=gfortran'
                              '--download-fblaslapack'
                              '--download-suitesparse'
                              '--with-debugging=1'
```

or ask the configure to download
```
python2 './configure' '--with-cc=gcc' '--with-cxx=g++'
                      '--with-fc=gfortran' '--download-fblaslapack'
                      '--download-suitesparse'
                      '--with-debugging=1' '--download-mpich'
```

To compile PETSc with MPI and HYPRE preconditioners
```
python2 './configure' --with-cc=gcc --with-cxx=g++ --with-fc=gfortran
--download-fblaslapack --download-suitesparse --download-hypre --download-mpich
```

After the configure, run the command
```
make PETSC_DIR=/home/rocha/source/petsc-3.X.X PETSC_ARCH=arch-linux2-c-debug all
```

## 2. PETSc tips

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
