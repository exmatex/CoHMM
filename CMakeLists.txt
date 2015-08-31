cmake_minimum_required(VERSION 2.6)

project(CoHMM_tba)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)

enable_language(C)
include(CMakeDetermineCCompiler)
enable_language(CXX)
include(CMakeDetermineCXXCompiler)

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

#for CoMD lib
include(CheckCCompilerFlag)
check_c_compiler_flag("-std=c99" HAVE_C99)
if(NOT HAVE_C99)
    message(FATAL_ERROR "flag '-std=c99' does not work")
endif(NOT HAVE_C99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -std=c99 -fPIC")

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" HAVE_CXX11)
if(NOT HAVE_CXX11)
    message(FATAL_ERROR "flag '-std=c++11' does not work")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -fPIC")


########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" OFF)

########################################################################
#Find external packages
########################################################################
find_package(HIREDIS REQUIRED)
include_directories(${HIREDIS_INCLUDE_DIRS})
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_package(SWIG)
if(SWIG_FOUND)
	INCLUDE(${SWIG_USE_FILE})
endif()

find_package(MPI)
if(MPI_CXX_FOUND)
	find_package(Libcircle)
	find_package(ChunksAndTasks)
endif()

find_package(TBB)
if(TBB_FOUND)
	find_package(CNC)
endif()

find_package(MKL)
find_package(CBLAS)
find_package(LAPACKE)
if(MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIR})
  set(LINALG_LIBRARIES ${MKL_LIBRARIES})
  add_definitions(-DHAVE_MKL)
elseif(CBLAS_FOUND AND LAPACKE_FOUND)
  include_directories(${CLBAS_INCLUDE_DIR})
  include_directories(${LAPACKE_INCLUDE_DIR})
  set(LINALG_LIBRARIES ${LAPACKE_LIBRARIES} ${CBLAS_LIBRARIES})
else()
  message(FATAL_ERROR "Neither mkl nor (blas + cblas + lapack + lapacke ) found!")
endif()

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFile)
foreach(HEADER math.h stdio.h stdlib.h string.h)
  check_include_file(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

include(CheckIncludeFileCXX)
foreach(HEADER algorithm cassert cfloat cmath cstdio cstdlib cstring ctime iostream list map vector)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
include(CheckLibraryExists)
foreach(FUNC sqrt)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)

######################################
# actual targets                     #
######################################
#CoMD
add_subdirectory(COMD_lib)

#RedisWrapper
add_subdirectory(RedisWrapper_lib)
#Everyone needs the RedisWrapper
include_directories(RedisWrapper_lib)

#CoHMM
add_subdirectory(COHMM_lib)
#Everyone needs CoHMM
include_directories(COHMM_lib)

#As targets are implemented, they leave the if(false)
if(false)
    #If SWIG found, build SWIG library
    if(SWIG_FOUND)
    	message(STATUS "SWIG Library Enabled")
    	set_source_files_properties(COHMM_lib/CoHMM_SwiftWrapper.i PROPERTIES CPLUSPLUS ON)
    	swig_add_module(cohmm_swiftt tcl8 COHMM_lib/CoHMM_SwiftWrapper.i ${COHMM_SOURCES})
    	target_link_libraries(cohmm_swiftt comd ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} -lrt)
    endif()

    #If Libcircle found, build Libcircle driver
    if(LIBCIRCLE_FOUND)
    	message(STATUS "Libcircle Driver Enabled")
    	add_subdirectory(2D_Libcircle)
    endif()

    if(MPI_CXX_FOUND)
    	message(STATUS "Naive MPI Driver Enabled")
        add_subdirectory(2D_MPI)
    endif()
endif()

#If Chunks And Tasks found, build Chunks and Tasks based driver
if(CHUNKSANDTASKS_FOUND)
    message(STATUS "Chunks and Tasks Driver Enabled")
    add_subdirectory(2D_CHT)
endif()

find_program(CHARMC NAMES charmc DOC "path to charmc")
if(CHARMC)
    message(STATUS "Charm++ Driver Enabled")
    add_subdirectory(2D_Charm++)
endif()

#If (Intel) CnC found, build CnC based driver
if(INTELCNC_FOUND)
    message(STATUS "Intel CnC Driver Enabled")
    add_subdirectory(2D_CnC)
endif()

if(OPENMP_FOUND)
	message(STATUS "OpenMP Driver Enabled")
    add_subdirectory(2D_OMP)
endif()
