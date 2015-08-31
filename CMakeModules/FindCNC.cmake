# - Find Intel's CnC
# Find the native Intel CnC headers and libraries.
#
#  INTELCNC_INCLUDE_DIRS - where to find hiredis.h, etc.
#  INTELCNC_LIBRARIES    - List of libraries when using hiredis.
#  INTELCNC_FOUND        - True if hiredis found.

# CnC needs TBB, so make sure that exists. The REQUIRED is forwarded
IF (NOT DEFINED TBB_FOUND)
    find_package(TBB)
ENDIF(NOT DEFINED TBB_FOUND)

# Look for the header file.
FIND_PATH(INTELCNC_INCLUDE_DIR NAMES cnc/cnc.h)

# Look for the library.
FIND_LIBRARY(INTELCNC_LIBRARY NAMES cnc libcnc )

# handle the QUIETLY and REQUIRED arguments and set INTELCNC_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(INTELCNC DEFAULT_MSG INTELCNC_LIBRARY INTELCNC_INCLUDE_DIR)

# Copy the results to the output variables.
IF(INTELCNC_FOUND)
  SET(INTELCNC_LIBRARIES ${INTELCNC_LIBRARY})
  SET(INTELCNC_INCLUDE_DIRS ${INTELCNC_INCLUDE_DIR})
ELSE(INTELCNC_FOUND)
  SET(INTELCNC_LIBRARIES)
  SET(INTELCNC_INCLUDE_DIRS)
ENDIF(INTELCNC_FOUND)

MARK_AS_ADVANCED(INTELCNC_INCLUDE_DIR INTELCNC_LIBRARY)

