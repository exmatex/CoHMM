cmake_minimum_required(VERSION 2.6)

include_directories(${TBB_INCLUDE_DIRS})
include_directories(${INTELCNC_INCLUDE_DIRS})

add_executable(2D_CnC 2D_CnC.cpp)
target_link_libraries(2D_CnC cohmm comd ${TBB_LIBRARIES} ${INTELCNC_LIBRARIES} ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} -lrt)

