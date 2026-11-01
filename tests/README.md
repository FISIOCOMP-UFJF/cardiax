# Cardiax test suite

```bash
./run_tests.sh                 # builds what changed and runs everything
./run_tests.sh --unit          # only subroutines and subproblems
./run_tests.sh --integration   # only complete problems
./run_tests.sh --no-gpu        # leaves out the AMGX tests
./run_tests.sh -R cardiac      # only the tests whose name matches
./run_tests.sh --help          # all the options
```

A failed test prints the assertion that failed and the end of the output of that case.
Each test file starts with a comment saying what it covers and why.

To run a single executable or a single case directly:

```bash
build/tests/bin/test_odes                  # all the cases of this file
build/tests/bin/test_odes --list           # names of the cases
build/tests/bin/test_odes "Land"           # the cases whose name contains "Land"
```

The output of each case stays in `build/tests/scratch/<test>/`, along with the files the
solvers wrote.

## Layout

| Path | Contents |
|------|----------|
| `tests/unit/` | Linked against the libraries. Subroutines and subproblems: cell models, element matrices, materials, circulation models. |
| `tests/integration/` | Run the executables in `build/app` on small problems, and check what they print and write. |
| `tests/support/` | `test_framework.hpp` (assertions), `test_main.cpp` (shared `main`), `fixtures.hpp` (generated meshes, running the solvers, reading HDF5). |
| `tests/golden/baseline.txt` | Regression baseline. See below. |

| Test | What it exercises |
|------|-------------------|
| `test_util` | command line, expressions, interpolation, input files, load control, parameters, time control, timers |
| `test_linalg` | PETSc vector, matrix and solver wrappers, 4th order tensors |
| `test_mesh` | every mesh format, writers, HDF5 reader, checkpoints, Purkinje mesh, graph colouring |
| `test_fem` | quadrature, shape functions, mappings, FE spaces, boundary values, L2 projection |
| `test_odes` | ODE integrators, all cell models, `Cells`, warm-up |
| `test_nls` | Newton and Newton with line search |
| `test_materials` | hyperelastic laws, active stress, stabilisation, mechanics utilities |
| `test_poisson_eikonal` | Poisson and Laplace assembly and convergence, eikonal |
| `test_mechanics_assembly` | tangent against finite differences of the residual, volumes, fibre stretch |
| `test_cardiac` | conductivity, stimuli, monodomain matrices, 0D circulation, circulation models against a lumped cavity, cell coupling, UVC transfer |
| `test_app_basic` | `poisson`, `l2projection`, `elasticity`, `eikonal` |
| `test_app_nonlinearelas` | `nonlinearelas`: closed forms, benchmarks, restart |
| `test_app_electrophysiology` | `monodomain` (also Purkinje), `bidomain` |
| `test_app_coupled` | `electromech`, `cardiac_cycle` with each circulation model, restart |
| `test_amgx` | the AMGX backend (label `gpu`), skipped if there is no GPU |

## Three kinds of expected result

**Known answers.** Exact integrals, manufactured solutions, closed forms, conservation laws.
A failure here means the code is wrong, or the test is.

**Regression baseline.** Where no closed form exists (an ionic model after a stimulus, the
state of a cardiac cycle), `CHECK_GOLDEN` compares against `tests/golden/baseline.txt`. These
values only say what the code produced when the baseline was taken: they detect changes and
do not prove the value is right. After an intended change in the results:

```bash
./run_tests.sh --update-golden
git diff tests/golden/baseline.txt   # review before keeping it
```

**Known bugs.** `TEST_CASE_KNOWN_BUG` documents a defect that exists today. The case states
the correct behaviour and is expected to fail; `run_tests.sh` lists these at the end, and they
do not turn the run red. When the code gets fixed the case starts passing, which is reported as
a failure until the marker is changed to `TEST_CASE`, and from then on it guards the fix.

## Adding a test

1. Start the file with a comment, in English, saying what is covered and why.
2. Write the cases with `TEST_CASE("...")`. For a solver input, generate a mesh with
   `write_tet_box` and friends, or copy a repository file with `copy_input`: solvers write
   next to their input, so they must not run on the original.
3. Add `cardiax_add_test(<file name> unit)` or `integration` to `tests/CMakeLists.txt`.

Things to keep in mind:

- `UpdatedLagrangian::solve()` and `pre_solve()` keep state in function-local statics. A process
  can only run one solve, so complete mechanics problems go through the executables.
- `petsc::Vector::operator=(double)` only accepts 0. With any other value it still zeroes the
  vector in a Release build, where the `assert` is compiled out. Use `set(value)`.
- Keep the problems small: the whole suite runs on every change.

