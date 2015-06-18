# - Find libcircle
# Find the native LIBCIRCLE headers and libraries.
#
#  LIBCIRCLE_INCLUDE_DIRS - where to find libcircle.h, etc.
#  LIBCIRCLE_LIBRARIES    - List of libraries when using libcircle.
#  LIBCIRCLE_FOUND        - True if libcircle found.

# Look for the header file.
FIND_PATH(LIBCIRCLE_INCLUDE_DIR NAMES libcircle.h PATH_SUFFIXES libcircle)

# Look for the library.
FIND_LIBRARY(LIBCIRCLE_LIBRARY NAMES circle )

# handle the QUIETLY and REQUIRED arguments and set LIBCIRCLE_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBCIRCLE DEFAULT_MSG LIBCIRCLE_LIBRARY LIBCIRCLE_INCLUDE_DIR)

# Copy the results to the output variables.
IF(LIBCIRCLE_FOUND)
  SET(LIBCIRCLE_LIBRARIES ${LIBCIRCLE_LIBRARY})
  SET(LIBCIRCLE_INCLUDE_DIRS ${LIBCIRCLE_INCLUDE_DIR})
ELSE(LIBCIRCLE_FOUND)
  SET(LIBCIRCLE_LIBRARIES)
  SET(LIBCIRCLE_INCLUDE_DIRS)
ENDIF(LIBCIRCLE_FOUND)

MARK_AS_ADVANCED(LIBCIRCLE_INCLUDE_DIR LIBCIRCLE_LIBRARY)
