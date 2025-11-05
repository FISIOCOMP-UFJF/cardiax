cmake_minimum_required(VERSION 2.8)
 
PROJECT(teste)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
#add_executable(convertFibersCell2Point convertFibersCell2Point.cpp)
add_executable(diffdispl diffdispl.cpp)
add_executable(vtu2msh vtu2msh.cpp)
add_executable(writeAPDvtu writeAPDvtu.cpp)
add_executable(writePartMeshVtu writePartMeshVtu.cpp)
add_executable(writeFibersVtu helper writeFibersVtu.cpp)
add_executable(refineHexMesh helper refineHexMesh.cpp)
 
message(${VTK_LIBRARIES})

#if(VTK_LIBRARIES)
#  target_link_libraries(convertFibersCell2Point ${VTK_LIBRARIES})
#endif()

#if(VTK_LIBRARIES)
#  target_link_libraries(diffdispl ${VTK_LIBRARIES})
#else()
#  target_link_libraries(diffdispl vtkHybrid)
#endif()

if(VTK_LIBRARIES)
  target_link_libraries(vtu2msh ${VTK_LIBRARIES})
else()
  target_link_libraries(vtu2msh vtkHybrid)
endif()

if(VTK_LIBRARIES)
  target_link_libraries(writeAPDvtu ${VTK_LIBRARIES})
else()
  target_link_libraries(writeAPDvtu vtkHybrid)
endif()

if(VTK_LIBRARIES)
  target_link_libraries(writePartMeshVtu ${VTK_LIBRARIES})
else()
  target_link_libraries(writePartMeshVtu vtkHybrid)
endif()

if(VTK_LIBRARIES)
  target_link_libraries(writeFibersVtu ${VTK_LIBRARIES})
else()
  target_link_libraries(writeFibersVtu vtkHybrid)
endif()

#if(VTK_LIBRARIES)
#  target_link_libraries(refineHexMesh ${VTK_LIBRARIES})
#else()
#  target_link_libraries(refineHexMesh vtkHybrid)
#endif()