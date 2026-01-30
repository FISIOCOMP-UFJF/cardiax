#!/bin/bash
set -e

# --- 1. ARGUMENT PROCESSING ---
FORCE_ALL=false
FORCE_HDF5=false
FORCE_ARMA=false
FORCE_AMGX=false
FORCE_PETSC=false
CUDA_ARCH=""

show_help() {
    echo "Usage: ./install_deps.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -f, --force-all       Forces recompilation of ALL dependencies."
    echo "  --force-hdf5          Forces recompilation of HDF5."
    echo "  --force-armadillo     Forces recompilation of Armadillo."
    echo "  --force-amgx          Forces recompilation of AMGX."
    echo "  --force-petsc         Forces recompilation of PETSc." 
    echo "  -h, --help            Shows this message. Hi! "
    echo ""
	
}

# Loop to read all passed flags
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--force-all)
            FORCE_ALL=true
            shift 
            ;;
        --force-hdf5)
            FORCE_HDF5=true
            shift
            ;;
        --force-armadillo)
            FORCE_ARMA=true
            shift
            ;;
        --force-amgx)
            FORCE_AMGX=true
            shift
            ;;
        --force-petsc)
            FORCE_PETSC=true
            shift
            ;;
	--cuda-arch)
            CUDA_ARCH="$2" #this is not working yet, but it will be necessary
	    shift 2
	    ;;

        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

echo "=========================================="
echo "   CARDIAX: DEPENDENCY INSTALLATION       "
if [ "$FORCE_ALL" = true ]; then
    echo "   (MODE: FORCE ALL)                      "
else
    echo "   (MODE: INCREMENTAL)            "
fi
echo "=========================================="

# 2. Initial Configuration
RECIPES_DIR="./recipes"

echo "[1/4] Checking Conda: "
eval "$(conda shell.bash hook)"

# 3. Setting up build environment
if ! conda info --envs | grep -q "cardiax_builder"; then
    echo "Creating build environment..."
    conda create -n cardiax_builder conda-build conda-index -c conda-forge -y
else
    echo "Environment 'cardiax_builder' OK."
fi

conda activate cardiax_builder

CONDA_BLD_PATH="$CONDA_PREFIX/conda-bld"
echo "Build path set to: $CONDA_BLD_PATH"


# 4. Smart Build Function
echo "[2/4] Verifying recipes:"

build_recipe_smart() {
    PKG_NAME=$1
    FORCE_THIS=$2 # Receives true/false if this specific package should be forced
    RECIPE_PATH="$RECIPES_DIR/$PKG_NAME"

    echo "--------------------------------------------------"
    echo "Analyzing: $PKG_NAME"

    # Ask Conda what the final file will be
    EXPECTED_OUTPUT=$(conda build "$RECIPE_PATH" \
        --output \
        -c local -c conda-forge -c nvidia \
        -m "$RECIPES_DIR/conda_build_config.yaml" \
        --output-folder "$CONDA_BLD_PATH" \
        2>/dev/null | tail -n 1 | tr -d '\r')

    # Decision Logic
    SHOULD_BUILD=false
    REASON=""

    if [ ! -f "$EXPECTED_OUTPUT" ]; then
        SHOULD_BUILD=true
        REASON="Package does not exist yet."
    elif [ "$FORCE_ALL" = true ]; then
        SHOULD_BUILD=true
        REASON="Flag --force-all activated."
    elif [ "$FORCE_THIS" = true ]; then
        SHOULD_BUILD=true
        REASON="Flag --force-$PKG_NAME activated."
    fi
	
    if [ "$SHOULD_BUILD" = true ]; then
        echo "[BUILDING] $PKG_NAME ($REASON)"
        
        # If forced, remove old file first to ensure clean build
        if [ -f "$EXPECTED_OUTPUT" ]; then
            rm "$EXPECTED_OUTPUT"
        fi

        # Light cleanup
        conda build purge > /dev/null 2>&1 || true

        conda build "$RECIPE_PATH" \
            -c local -c conda-forge -c nvidia \
            -m "$RECIPES_DIR/conda_build_config.yaml" \
            --output-folder "$CONDA_BLD_PATH"
        
        echo "   -> Updating local index..."
        python -m conda_index "$CONDA_BLD_PATH" > /dev/null 2>&1
    else
        echo "[SKIPPING] $PKG_NAME already exists."
    fi
}

build_recipe_smart "hdf5_custom"      "$FORCE_HDF5"
build_recipe_smart "armadillo_custom" "$FORCE_ARMA"

if [ -n "$CUDA_ARCH" ]; then
    echo "CUDA architectures set to: $CUDA_ARCH"
    export CARDIAX_CUDA_ARCH="$CUDA_ARCH"
fi

build_recipe_smart "amgx_custom"      "$FORCE_AMGX"
build_recipe_smart "petsc_custom"     "$FORCE_PETSC"

# 6. Final Indexing
echo "[3/4] Finalizing index of dependencies:"
python -m conda_index "$CONDA_BLD_PATH"

# 7. Create Final Environment
echo "[4/4] Recreating final environment 'cardiax_env':"
conda deactivate

conda remove -n cardiax_env --all -y > /dev/null 2>&1 || true

conda create -n cardiax_env \
    -c "file://$CONDA_BLD_PATH" \
    -c conda-forge -c nvidia \
    hdf5_custom armadillo_custom amgx_custom petsc_custom \
    cuda-nvtx cuda-nvtx-dev cuda-libraries-dev cuda-cudart-dev "cuda-version=12" \
    cmake make gxx_linux-64 gcc_linux-64 gfortran_linux-64 \
    python=3.10 mpich pkg-config -y

echo "=========================================="
echo "   SUCCESS! All done. Run the ./build.sh  "
echo "           to compile Cardiax             "
echo "=========================================="
