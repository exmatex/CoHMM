cmake_minimum_required(VERSION 2.6)

add_custom_command(
  OUTPUT CoMD_info.h
  COMMAND "${CMAKE_SOURCE_DIR}/COMD_lib/generate_info_header" "${CMAKE_C_COMPILER}"
  "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}"
  "${CMAKE_EXE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}"
  DEPENDS  "${CMAKE_SOURCE_DIR}/COMD_lib/generate_info_header"
  COMMENT "Generating CoMD_info.h ...")
add_custom_target(comd_header DEPENDS CoMD_info.h)

include_directories(${CMAKE_BINARY_DIR}/COMD_lib)
file(GLOB COMD_SOURCES *.c)
add_library(comd ${COMD_SOURCES})
add_dependencies(comd comd_header)
