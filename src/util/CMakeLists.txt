# CmakeLists.txt in src/util dir
# Adds a library called Util from the source files

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

target_include_directories(Util PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)


target_link_libraries(Util PUBLIC
    ext_libs
)


add_executable(testParameters test/test_parameters.cpp)

target_link_libraries(testParameters PRIVATE Util)