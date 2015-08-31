# - Find hiredis
# Find the native HIREDIS headers and libraries.
#
#  HIREDIS_INCLUDE_DIRS - where to find hiredis.h, etc.
#  HIREDIS_LIBRARIES    - List of libraries when using hiredis.
#  HIREDIS_FOUND        - True if hiredis found.

# Look for the header file.
FIND_PATH(HIREDIS_INCLUDE_DIR NAMES hiredis.h PATH_SUFFIXES hiredis)

# Look for the library.
FIND_LIBRARY(HIREDIS_LIBRARY NAMES hiredis libhiredis )

# handle the QUIETLY and REQUIRED arguments and set HIREDIS_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HIREDIS DEFAULT_MSG HIREDIS_LIBRARY HIREDIS_INCLUDE_DIR)

# Copy the results to the output variables.
IF(HIREDIS_FOUND)
  SET(HIREDIS_LIBRARIES ${HIREDIS_LIBRARY})
  SET(HIREDIS_INCLUDE_DIRS ${HIREDIS_INCLUDE_DIR})
ELSE(HIREDIS_FOUND)
  SET(HIREDIS_LIBRARIES)
  SET(HIREDIS_INCLUDE_DIRS)
ENDIF(HIREDIS_FOUND)

MARK_AS_ADVANCED(HIREDIS_INCLUDE_DIR HIREDIS_LIBRARY)
