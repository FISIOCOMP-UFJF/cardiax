# CmakeLists.txt in src dir
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

add_executable(testHDF5 ../mesh/test/test_hdf5.cpp)
target_link_libraries(testHDF5 Mesh FEM Util LinAlg ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES}  ${PETSC_LIBRARIES})

add_executable(testMesh ../mesh/test/test_mesh.cpp)
target_link_libraries(testMesh Mesh FEM Util LinAlg ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES}  ${PETSC_LIBRARIES})

add_executable(testXML ../mesh/test/test_xml.cpp)
target_link_libraries(testXML Mesh FEM Util LinAlg ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES}  ${PETSC_LIBRARIES})

add_executable(testVTK ../mesh/test/test_vtk.cpp)
target_link_libraries(testVTK Mesh FEM Util LinAlg ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES}  ${PETSC_LIBRARIES})

