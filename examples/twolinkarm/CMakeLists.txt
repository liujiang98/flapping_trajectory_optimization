cmake_minimum_required (VERSION 3.10)
project(twolinkarm LANGUAGES CXX)
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} )

SET(CUSTOM_RBDL_PATH "" CACHE PATH "Path to specific RBDL Installation")
FIND_PACKAGE (RBDL COMPONENTS URDFREADER REQUIRED)

add_executable(${PROJECT_NAME} $<$<BOOL:$<TARGET_EXISTS:PSOPT_SNOPT_interface>>:$<TARGET_OBJECTS:PSOPT_SNOPT_interface>> twolinkarm.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt adolc Eigen3::Eigen ${RBDL_LIBRARY} ${RBDL_URDFREADER_LIBRARY})

if(${WITH_SNOPT_INTERFACE})
	enable_language(Fortran)
	target_link_libraries(${PROJECT_NAME} PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)
endif()
