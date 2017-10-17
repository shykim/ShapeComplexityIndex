cmake_minimum_required (VERSION 2.8)
project (Cal_Complexity)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
	INCLUDE(${GenerateCLP_USE_FILE})
ENDIF(GenerateCLP_FOUND)
generateCLP(CLP_XML ShapeComplexityIndex.xml)

#find_path(EIGEN3 Dense)
#include_directories( ${EIGEN3} )

# MESH library
add_library(Mesh SHARED IMPORTED)

IF(APPLE)
	set_target_properties(Mesh PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libMesh.1.0.dylib)
ELSE()
	set_target_properties(Mesh PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libMesh.so.1.0)
ENDIF(APPLE)	
	
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#CONFIGURE_FILE(ShapeComplexityIndex.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Cal_Complex_opt.x.cpp)

CONFIGURE_FILE(surf_reg_model_left.obj template.obj COPYONLY)

# executable file
add_executable(Cal_Complex Cal_Complex_opt.x.cpp Calculators.h Calculators.cpp transportSimplex.h ${CLP_XML})
target_link_libraries(Cal_Complex Mesh)


