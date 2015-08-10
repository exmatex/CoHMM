cmake_minimum_required(VERSION 2.6)

include_directories(${TBB_INCLUDE_DIRS})
include_directories(${INTELCNC_INCLUDE_DIRS})

add_definitions(-DCNC_DIST)

add_executable(2D_CnC CoHMM_CnC.cpp 2DKriging.cpp kriging.cpp redisBuckets.cpp cncShmem.cpp CoHMM_CnCaD.cpp)
target_link_libraries(2D_CnC comd ${TBB_LIBRARIES} ${INTELCNC_LIBRARIES} ${LINALG_LIBRARIES} ${HIREDIS_LIBRARIES} -lrt)
