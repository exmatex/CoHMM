cmake_minimum_required(VERSION 2.6)

add_library(RedisWrapper RedisWrapper.cpp)
target_link_libraries(RedisWrapper ${HIREDIS_LIBRARIES})
