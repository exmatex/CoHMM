cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_COMPILER "${CHARMC}")
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_executable(2D_Charm++ 2D_Charm++.cpp)

file(GLOB COHMM_DAD_CHARM_MODULES *.ci)
foreach(_CI ${COHMM_DAD_CHARM_MODULES})
	get_filename_component(_TRUNC "${_CI}" NAME_WE)
	add_custom_command(OUTPUT ${_TRUNC}.decl.h ${_TRUNC}.def.h
		COMMAND ${CHARMC} ${_CI}
		COMMENT "Compiling ${_CI}")
	add_custom_target(${_TRUNC} DEPENDS
		${_TRUNC}.decl.h
		${_TRUNC}.def.h)
	add_dependencies(2D_Charm++ ${_TRUNC})
endforeach()

target_link_libraries(2D_Charm++ cohmm comd ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} RedisWrapper  -lrt)
