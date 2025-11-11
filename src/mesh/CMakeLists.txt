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

# [NOVO] Exporta os headers do Mesh
target_include_directories(Mesh PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# [NOVO] Declara as dependências da biblioteca Mesh
# Baseado nos seus arquivos (writer_hdf5.cpp) e testes,
# Mesh depende de Util, LinAlg e libs externas.
target_link_libraries(Mesh PUBLIC
    Util
    LinAlg
    ext_libs # Nosso alvo que tem HDF5, PETSc, Armadillo
)
# NOTA: Você tem "writer_vtk.cpp". Você deve adicionar
# find_package(VTK REQUIRED) no seu CMakeLists.txt principal
# e adicionar ${VTK_LIBRARIES} e ${VTK_INCLUDE_DIRS} ao 'ext_libs'.


# --- Binaries for tests ---

# [NOVO] Para evitar repetição, criamos um alvo de interface
# que agrupa as dependências comuns de todos os testes neste arquivo.
add_library(mesh_test_deps INTERFACE)
target_link_libraries(mesh_test_deps INTERFACE
    Mesh # Dependência principal
    FEM  # Os testes também precisam de FEM
)
# Como Mesh e FEM já linkam 'ext_libs', 'Util', e 'LinAlg'
# como PUBLIC, 'mesh_test_deps' agora contém TUDO
# que os testes precisam (PETSc, HDF5, Armadillo, etc.)


# Adiciona os executáveis de teste
add_executable(testHDF5 ../mesh/test/test_hdf5.cpp)
add_executable(testMesh ../mesh/test/test_mesh.cpp)
add_executable(testXML ../mesh/test/test_xml.cpp)
add_executable(testVTK ../mesh/test/test_vtk.cpp)

# [MODIFICADO] Linka todos os testes contra o alvo de interface
# Em vez de 4 chamadas longas e repetidas, temos 4 chamadas curtas.
target_link_libraries(testHDF5 PRIVATE mesh_test_deps)
target_link_libraries(testMesh PRIVATE mesh_test_deps)
target_link_libraries(testXML PRIVATE mesh_test_deps)
target_link_libraries(testVTK PRIVATE mesh_test_deps)

# --- REMOVIDO ---
# As 4 chamadas target_link_libraries(...) longas e repetitivas
# foram substituídas pela abordagem acima.