cmake_minimum_required(VERSION 2.6)

file(GLOB COHMM_SOURCES *.cpp)
include_directories(${CMAKE_SOURCE_DIR}/COMD_lib)

add_library(cohmm ${COHMM_SOURCES})
target_link_libraries(cohmm comd RedisWrapper ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} -lrt)

if(false)
    #If SWIG found, build SWIG library
    if(SWIG_FOUND)
    	message(STATUS "SWIG Library Enabled")
    	set_source_files_properties(CoHMM_SwiftWrapper.i PROPERTIES CPLUSPLUS ON)
    	swig_add_module(cohmm_swiftt tcl8 CoHMM_SwiftWrapper.i ${COHMM_SOURCES})
    	target_link_libraries(cohmm_swiftt comd ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} -lrt)
    endif()
endif()
