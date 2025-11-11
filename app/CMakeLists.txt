# Simple CMakeLists for building the applications
# Bernardo M. Rocha

# --- REMOVIDO ---
# set(CORELIBS ...)
# set(BOOSTLIBS ...)
# A linkagem agora é gerenciada por alvos e transitividade.
# As bibliotecas externas (BLAS, HDF5, MPI, Boost, PETSc)
# estão todas agrupadas no alvo 'ext_libs'.

set(PROGRAMS
    l2projection
    poisson
    elasticity
    nonlinearelas
    monodomain
    bidomain
    electromech
)

# Build them all
foreach(program ${PROGRAMS})
    # [SEM MUDANÇA] Define o executável
    add_executable(${program} ${program}_solver.cpp)

    # [MODIFICADO] Linkagem moderna
    # Nós linkamos contra os alvos de biblioteca de "nível superior"
    # do nosso projeto (provavelmente 'Cardiac' e 'PDE').
    #
    # Como 'Cardiac' e 'PDE' foram linkados (em seus CMakeLists)
    # como PUBLIC às suas dependências (FEM, ODE, LinAlg, Util)
    # e ao 'ext_libs' (que tem PETSc, HDF5, Boost, MPI, BLAS, etc.),
    # o CMake lida com TODAS as dependências transitivas.
    #
    # O executável não precisa saber sobre PETSc, HDF5, ou
    # mesmo 'Util'. Ele só precisa saber sobre 'Cardiac' e 'PDE'.
    target_link_libraries(${program} PRIVATE
        Cardiac
        PDE
    )

    # --- REMOVIDO ---
    # A linha antiga era:
    # link_directories(${Boost_LIBRARIES})
    # target_link_libraries(${program} ${CORELIBS} ${Boost_LIBRARIES} ${BOOSTLIBS} ${PETSC_LIBRARIES})
endforeach(program)

# --- NOTA SOBRE PARDISO ---
# Seu comentário sobre PARDISO (pardiso500-GNU481-X86-64)
# deve ser gerenciado no CMakeLists.txt principal,
# encontrando-o com find_library() e adicionando-o
# ao alvo 'ext_libs'. Não o coloque "hard-coded" aqui.