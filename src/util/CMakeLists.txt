# CmakeLists.txt in src/util dir
# Adds a library called Util from the source files

# [MODIFICADO] Lista de fontes da biblioteca (SEM o teste)
add_library(Util
    command_line_args.cpp
    interp.cpp
    input_file.cpp
    log.cpp
    load_control.cpp
    parameter.cpp
    parameters.cpp
    time_parameters.cpp
    pugixml.cpp
    timer.cpp
    util.cpp
)

target_include_directories(Util PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# [NOVO] Adiciona as dependências do Util
# Util é a sua biblioteca base, ela precisa das bibliotecas externas
# (PETSc, Boost, etc.) que estão no nosso alvo 'ext_libs'.
# Usamos PUBLIC para que qualquer alvo que linkar 'Util'
# (como 'FEM', 'ODE', etc.) receba 'ext_libs' transitivamente.
target_link_libraries(Util PUBLIC
    ext_libs
)


# --- Testes ---
# [MODIFICADO] Define o executável de teste
# O caminho está corrigido (assumindo que 'test' é um subdiretório)
add_executable(testParameters test/test_parameters.cpp)

# [MODIFICADO] Linka o teste
# O teste depende da biblioteca 'Util'.
# Como 'Util' linka 'ext_libs' como PUBLIC, 'testParameters'
# receberá todas as bibliotecas (PETSc, Armadillo, etc.)
# automaticamente, por transitividade.
target_link_libraries(testParameters PRIVATE Util)