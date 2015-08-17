cmake_minimum_required(VERSION 2.6)

add_executable(2D_OMP 2D_OMP.cpp)
target_link_libraries(2D_OMP cohmm comd ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} RedisWrapper -lrt)
