
cmake_minimum_required(VERSION 2.8)

PROJECT(SupernovaSeg)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(SupernovaSeg MACOSX_BUNDLE SupernovaSeg _database _Supernova)

if(VTK_LIBRARIES)
  target_link_libraries(SupernovaSeg ${VTK_LIBRARIES})
else()
  target_link_libraries(SupernovaSeg vtkHybrid vtkWidgets)
endif()

