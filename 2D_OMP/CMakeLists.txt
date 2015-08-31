cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(2D_OMP CoHMM_OMP.cpp 2D_OMP.cpp)
target_link_libraries(2D_OMP comd cohmm RedisWrapper ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} -lrt)
