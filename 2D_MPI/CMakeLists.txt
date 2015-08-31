cmake_minimum_required(VERSION 2.6)

include_directories(${MPI_CXX_INCLUDE_PATH})

add_executable(2D_MPI 2D_MPI.cpp)

if(MPI_COMPILE_FLAGS)
	set_target_properties(2D_MPI PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
	set_target_properties(2D_MPI PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()


target_link_libraries(2D_MPI cohmm comd ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} ${MPI_CXX_LIBRARIES} RedisWrapper -lrt)
