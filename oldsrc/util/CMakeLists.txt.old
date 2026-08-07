# CmakeLists.txt in src dir
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
    test/test_parameters.cpp)

# Programs
#if(BUILD_TESTS)
#    add_executable(testInterpolation test_interp.cpp)
#    target_link_libraries(testInterpolation m z Util LinAlg)
#endif(BUILD_TESTS)

add_executable(testParameters ../util/test/test_parameters.cpp)
target_link_libraries(testParameters Util)