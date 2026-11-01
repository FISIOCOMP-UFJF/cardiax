# Cardiax test suite. See tests/README.md.
#
#   unit/         subroutines and subproblems, linked against the libraries
#   integration/  complete problems, running the executables in build/app
#
# Everything here is EXCLUDE_FROM_ALL (see the top-level CMakeLists.txt), so
# ./build.sh does not pay for it. ./run_tests.sh builds the `cardiax_tests` target.

add_library(cardiax_test_support STATIC
  support/test_main.cpp
  support/fixtures.cpp
)
target_include_directories(cardiax_test_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/support)
target_compile_definitions(cardiax_test_support PUBLIC
  CARDIAX_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
  CARDIAX_APP_DIR="${PROJECT_BINARY_DIR}/app"
  CARDIAX_SCRATCH_ROOT="${PROJECT_BINARY_DIR}/tests/scratch"
)
target_link_libraries(cardiax_test_support PUBLIC Cardiac PDE kdtree_uvc)

add_custom_target(cardiax_tests)

# cardiax_add_test(<name> <unit|integration> [TIMEOUT s] [ARGS ...] [LABELS ...])
# The source is <kind>/<name>.cpp. Two threads per process: results do not depend
# on it, and ctest runs several tests at once.
function(cardiax_add_test name kind)
  cmake_parse_arguments(T "" "TIMEOUT;SOURCE" "ARGS;LABELS" ${ARGN})
  if(NOT T_TIMEOUT)
    set(T_TIMEOUT 300)
  endif()
  if(NOT T_SOURCE)
    set(T_SOURCE ${name})
  endif()

  if(NOT TARGET ${T_SOURCE}_exe)
    add_executable(${T_SOURCE}_exe ${kind}/${T_SOURCE}.cpp)
    set_target_properties(${T_SOURCE}_exe PROPERTIES
      OUTPUT_NAME ${T_SOURCE}
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
    target_link_libraries(${T_SOURCE}_exe PRIVATE cardiax_test_support)
    add_dependencies(cardiax_tests ${T_SOURCE}_exe)
  endif()

  add_test(NAME ${name} COMMAND ${T_SOURCE}_exe ${T_ARGS})
  set_tests_properties(${name} PROPERTIES
    LABELS "${kind};${T_LABELS}"
    TIMEOUT ${T_TIMEOUT}
    ENVIRONMENT "OMP_NUM_THREADS=2")
endfunction()

# The integration tests run these
add_dependencies(cardiax_tests
  l2projection poisson elasticity nonlinearelas monodomain bidomain
  cardiac_cycle electromech eikonal)

# --- Unit tests ---
cardiax_add_test(test_util unit)
cardiax_add_test(test_linalg unit)
cardiax_add_test(test_mesh unit)
cardiax_add_test(test_fem unit)
cardiax_add_test(test_odes unit TIMEOUT 600)
cardiax_add_test(test_nls unit)
cardiax_add_test(test_materials unit)
cardiax_add_test(test_poisson_eikonal unit)
cardiax_add_test(test_cardiac unit)
cardiax_add_test(test_mechanics_assembly unit)

# --- Integration tests: complete problems, through the executables ---
cardiax_add_test(test_app_basic integration)
cardiax_add_test(test_app_nonlinearelas integration TIMEOUT 600)
cardiax_add_test(test_app_electrophysiology integration TIMEOUT 900)
cardiax_add_test(test_app_coupled integration TIMEOUT 1800)

# --- GPU tests: the AMGX backend. Skipped at run time if there is no GPU. ---
cardiax_add_test(test_amgx integration LABELS gpu)
set_tests_properties(test_amgx PROPERTIES RESOURCE_LOCK gpu)
