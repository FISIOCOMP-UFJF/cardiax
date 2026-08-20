# CmakeLists.txt in src/mesh dir
# Adds a library called Mesh from the source files

add_library(Mesh
    mesh.cpp
    box_mesh.cpp
    purkinje_mesh.cpp
    writer.cpp
    writer_vtk.cpp
    writer_hdf5.cpp
    gmshio.cpp
    graph_coloring.cpp
)

target_include_directories(Mesh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Mesh PUBLIC
    Util
    LinAlg
    ext_libs # Nosso alvo que tem HDF5, PETSc, Armadillo
)


# --- Binaries for tests ---
add_library(mesh_test_deps INTERFACE)
target_link_libraries(mesh_test_deps INTERFACE
    Mesh 
    FEM  
)


add_executable(testHDF5 ../mesh/test/test_hdf5.cpp)
add_executable(testMesh ../mesh/test/test_mesh.cpp)
add_executable(testXML ../mesh/test/test_xml.cpp)
add_executable(testVTK ../mesh/test/test_vtk.cpp)

target_link_libraries(testHDF5 PRIVATE mesh_test_deps)
target_link_libraries(testMesh PRIVATE mesh_test_deps)
target_link_libraries(testXML PRIVATE mesh_test_deps)
target_link_libraries(testVTK PRIVATE mesh_test_deps)
