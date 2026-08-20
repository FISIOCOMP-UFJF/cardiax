add_library(NLS
    nonlinear_solver.cpp
    bfgs.cpp
    newton.cpp
    newton_ls.cpp
)

target_include_directories(NLS PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)


target_link_libraries(NLS PUBLIC
    LinAlg
    ext_libs # Fornece PETSc, AMGX, gomp, zlib, etc.
)

add_executable(testSimple EXCLUDE_FROM_ALL
    test/simple.cpp
    test/test_simple.cpp
)

target_link_libraries(testSimple PRIVATE
    NLS
    ODE
    Util
)
