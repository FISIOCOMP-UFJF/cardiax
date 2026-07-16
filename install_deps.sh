#!/bin/bash
set -e

# --- 1. ARGUMENT PROCESSING ---
FORCE_ALL=false
FORCE_HDF5=false
FORCE_ARMA=false
FORCE_AMGX=false
FORCE_PETSC=false
CUDA_ARCH=""

IGNORE_HDF5=false
IGNORE_ARMA=false
IGNORE_AMGX=false
IGNORE_PETSC=false

show_help() {
    echo "Usage: ./install_deps.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -f, --force-all       Forces recompilation of ALL dependencies."
    echo "  --force-hdf5          Forces recompilation of HDF5."
    echo "  --force-armadillo     Forces recompilation of Armadillo."
    echo "  --force-amgx          Forces recompilation of AMGX."
    echo "  --force-petsc         Forces recompilation of PETSc." 
    echo "  --ignore-hdf5         Skips HDF5 build and installation."
    echo "  --ignore-armadillo    Skips Armadillo build and installation."
    echo "  --ignore-amgx         Skips AMGX build and installation."
    echo "  --ignore-petsc        Skips PETSc build and installation."
    echo "  --cuda-arch           Specifies the CUDA architecture."
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
        --ignore-hdf5)
            IGNORE_HDF5=true
            shift
            ;;
        --ignore-armadillo)
            IGNORE_ARMA=true
            shift
            ;;
        --ignore-amgx)
            IGNORE_AMGX=true
            shift
            ;;
        --ignore-petsc)
            IGNORE_PETSC=true
            shift
            ;;
        --cuda-arch)
            CUDA_ARCH="$2"
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
    echo "   (MODE: INCREMENTAL FAST)               "
fi
echo "=========================================="

# 2. Initial Configuration
RECIPES_DIR="./recipes"
MARKERS_DIR="./.build_markers"

# Cria a pasta de marcadores se ela não existir
mkdir -p "$MARKERS_DIR"

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


# 4. Smart Build Function (Otimizada com Marcadores)
echo "[2/4] Verifying recipes:"

build_recipe_smart() {
    PKG_NAME=$1
    FORCE_THIS=$2 # Receives true/false if this specific package should be forced
    RECIPE_PATH="$RECIPES_DIR/$PKG_NAME"
    MARKER_FILE="$MARKERS_DIR/${PKG_NAME}.done"

    echo "--------------------------------------------------"
    echo "Analyzing: $PKG_NAME"

    # Decision Logic baseada em arquivos de texto rápidos
    SHOULD_BUILD=false
    REASON=""

    if [ ! -f "$MARKER_FILE" ]; then
        SHOULD_BUILD=true
        REASON="Marker not found (Not built yet)."
    elif [ "$FORCE_ALL" = true ]; then
        SHOULD_BUILD=true
        REASON="Flag --force-all activated."
    elif [ "$FORCE_THIS" = true ]; then
        SHOULD_BUILD=true
        REASON="Flag --force-$PKG_NAME activated."
    fi
    
    if [ "$SHOULD_BUILD" = true ]; then
        echo "[BUILDING] $PKG_NAME ($REASON)"
        
        # Remove o marcador antigo caso seja uma recompilação forçada
        if [ -f "$MARKER_FILE" ]; then
            rm -f "$MARKER_FILE"
        fi

        # Light cleanup
        conda build purge > /dev/null 2>&1 || true

        conda build "$RECIPE_PATH" \
            -c local -c conda-forge -c nvidia \
            -m "$RECIPES_DIR/conda_build_config.yaml" \
            --output-folder "$CONDA_BLD_PATH"
        
        echo "   -> Updating local index..."
        python -m conda_index "$CONDA_BLD_PATH" > /dev/null 2>&1

        # Cria o arquivo marcador para pular essa compilação no futuro.
        touch "$MARKER_FILE"
        echo "[SUCCESS] Marker created for $PKG_NAME"
    else
        echo "[SKIPPING] $PKG_NAME already exists (Fast check via marker)."
    fi
}

# --- Execução Condicional dos Builds ---

if [ "$IGNORE_HDF5" = false ]; then
    build_recipe_smart "hdf5_custom" "$FORCE_HDF5"
else
    echo "[IGNORED] hdf5_custom"
fi

if [ "$IGNORE_ARMA" = false ]; then
    build_recipe_smart "armadillo_custom" "$FORCE_ARMA"
else
    echo "[IGNORED] armadillo_custom"
fi

if [ -n "$CUDA_ARCH" ]; then
    echo "CUDA architectures set to: $CUDA_ARCH"
    export CARDIAX_CUDA_ARCH="$CUDA_ARCH"
fi

if [ "$IGNORE_AMGX" = false ]; then
    build_recipe_smart "amgx_custom" "$FORCE_AMGX"
else
    echo "[IGNORED] amgx_custom"
fi

if [ "$IGNORE_PETSC" = false ]; then
    build_recipe_smart "petsc_custom" "$FORCE_PETSC"
else
    echo "[IGNORED] petsc_custom"
fi

# 6. Final Indexing
echo "[3/4] Finalizing index of dependencies:"
python -m conda_index "$CONDA_BLD_PATH"

# 7. Create Final Environment
echo "[4/4] Recreating final environment 'cardiax_env':"
conda deactivate

conda remove -n cardiax_env --all -y > /dev/null 2>&1 || true

# Constrói a string de dependências apenas com o que não foi ignorado
DEPS=""
if [ "$IGNORE_HDF5" = false ]; then DEPS="$DEPS hdf5_custom"; fi
if [ "$IGNORE_ARMA" = false ]; then DEPS="$DEPS armadillo_custom"; fi
if [ "$IGNORE_AMGX" = false ]; then DEPS="$DEPS amgx_custom"; fi
if [ "$IGNORE_PETSC" = false ]; then DEPS="$DEPS petsc_custom"; fi

conda create -n cardiax_env \
    -c "file://$CONDA_BLD_PATH" \
    -c conda-forge -c nvidia \
    $DEPS \
    cuda-nvtx cuda-nvtx-dev cuda-libraries-dev cuda-cudart-dev "cuda-version=12" \
    cmake make gxx_linux-64 gcc_linux-64 gfortran_linux-64 \
    python=3.10 mpich pkg-config -y

echo "=========================================="
echo "   SUCCESS! All done. Run the ./build.sh  "
echo "           to compile Cardiax             "
echo "=========================================="