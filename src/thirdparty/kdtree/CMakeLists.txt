add_library(kdtree STATIC kdtree.cpp)

# anyone linking kdtree gets this include dir automatically
target_include_directories(kdtree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# needed if kdtree ends up inside a shared library (libcardiax.so, etc.)
set_target_properties(kdtree PROPERTIES POSITION_INDEPENDENT_CODE ON)