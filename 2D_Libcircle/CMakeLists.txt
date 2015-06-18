cmake_minimum_required(VERSION 2.6)

include_directories(${LIBCIRCLE_INCLUDE_DIRS})
add_executable(2D_Libcircle 2D_Libcircle.cpp )
target_link_libraries(2D_Libcircle cohmm comd ${LIBCIRCLE_LIBRARY} ${MPI_LIBRARIES} ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} -lrt)
if(MPI_COMPILE_FLAGS)
	set_target_properties(2D_Libcircle PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
	set_target_properties(2D_Libcircle PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

