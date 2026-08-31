# Simple CMakeLists for building the applications
# Bernardo M. Rocha

set(PROGRAMS
    l2projection
    poisson
    elasticity
    nonlinearelas
    monodomain
    bidomain
    electromech
    coupled_electromech
)

# Build them all
foreach(program ${PROGRAMS})
    add_executable(${program} ${program}_solver.cpp)

    target_link_libraries(${program} PRIVATE
        Cardiac
        PDE
    )

endforeach(program)
