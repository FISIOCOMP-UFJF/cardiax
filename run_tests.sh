#!/bin/bash
# Builds and runs the Cardiax test suite. See tests/README.md.

usage() {
  cat <<EOF
Usage: ./run_tests.sh [OPTIONS]

With no options, builds everything and runs ALL the tests.

Options:
  --unit             only the unit tests (subroutines and subproblems)
  --integration      only the integration tests (complete problems)
  --no-gpu           leave out the tests that use AMGX on the GPU
  -R, --regex <re>   only the tests whose name matches <re>
  -j <n>             tests running at the same time (default: 3)
  --build-jobs <n>   compilation jobs (default: all cores)
  --no-build         do not build, run what is already built
  --update-golden    rewrite tests/golden/baseline.txt with the current results
                     (only after an INTENDED change; review the diff)
  -V, --verbose      show the output of the tests that pass as well
  -h, --help         this message
EOF
}

LABEL=""
EXCLUDE=""
REGEX=""
JOBS=3
BUILD_JOBS=$(nproc)
DO_BUILD=1
GOLDEN=0
VERBOSE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --unit)          LABEL="unit"; shift ;;
    --integration)   LABEL="integration"; shift ;;
    --no-gpu)        EXCLUDE="gpu"; shift ;;
    -R|--regex)      REGEX="$2"; shift 2 ;;
    -j)              JOBS="$2"; shift 2 ;;
    --build-jobs)    BUILD_JOBS="$2"; shift 2 ;;
    --no-build)      DO_BUILD=0; shift ;;
    --update-golden) GOLDEN=1; shift ;;
    -V|--verbose)    VERBOSE="-V"; shift ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 2 ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

# The same environment as build.sh
eval "$(conda shell.bash hook)" || exit 1
conda activate cardiax_env || { echo "Run ./install_deps.sh first."; exit 1; }

# Not called BUILD: activating the conda compilers exports a variable with that name
if [[ $DO_BUILD -eq 1 ]]; then
  if [[ ! -f build/Makefile ]]; then
    ./build.sh || exit 1
  fi
  make -C build -j"$BUILD_JOBS" cardiax_tests || exit 1
fi

CTEST_ARGS=(--test-dir build --output-on-failure -j "$JOBS")
[[ -n "$LABEL" ]]   && CTEST_ARGS+=(-L "$LABEL")
[[ -n "$EXCLUDE" ]] && CTEST_ARGS+=(-LE "$EXCLUDE")
[[ -n "$REGEX" ]]   && CTEST_ARGS+=(-R "$REGEX")
[[ -n "$VERBOSE" ]] && CTEST_ARGS+=("$VERBOSE")

if [[ $GOLDEN -eq 1 ]]; then
  export CARDIAX_UPDATE_GOLDEN=1
  # The tests rewrite one shared file
  CTEST_ARGS+=(-j 1)
fi

ctest "${CTEST_ARGS[@]}"
STATUS=$?

# Known bugs are not failures, but they must stay visible
LOG="build/Testing/Temporary/LastTest.log"
if [[ -f "$LOG" ]] && grep -q "KNOWN BUG" "$LOG"; then
  echo
  echo "Known bugs in the code (documented by the tests, not counted as failures):"
  grep -A1 "KNOWN BUG" "$LOG" | grep -v '^--' | sed 's/^/  /'
fi

exit $STATUS
