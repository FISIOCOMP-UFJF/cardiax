option(UVC_BUILD_TESTS "Build teste_kdtree_uvc" ON)

add_library(kdtree_uvc STATIC
  kdtree_uvc.cpp)

target_include_directories(kdtree_uvc PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(kdtree_uvc PUBLIC
  ReaderHDF5
  kdtree)

if(UVC_BUILD_TESTS)
  add_executable(teste_kdtree_uvc teste_kdtree_uvc.cpp)
  target_link_libraries(teste_kdtree_uvc PRIVATE kdtree_uvc)
endif()
