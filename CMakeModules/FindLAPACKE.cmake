# - Find lapacke
# Find the native LAPACKE headers and libraries.
#
#  LAPACKE_INCLUDE_DIRS - where to find lapacke.h, etc.
#  LAPACKE_LIBRARIES    - List of libraries when using lapacke.
#  LAPACKE_FOUND        - True if lapacke found.

# Look for the header file.
FIND_PATH(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATH_SUFFIXES lapacke)

# Look for the library.
FIND_LIBRARY(LAPACKE_LIBRARY NAMES lapacke reflapacke)

# lapacke needs lapack
FIND_PACKAGE(LAPACK)

# handle the QUIETLY and REQUIRED arguments and set LAPACKE_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LAPACKE DEFAULT_MSG LAPACKE_LIBRARY LAPACK_LIBRARIES LAPACKE_INCLUDE_DIR)

# Copy the results to the output variables.
IF(LAPACKE_FOUND)
  SET(LAPACKE_LIBRARIES ${LAPACKE_LIBRARY} ${LAPACK_LIBRARIES})
  SET(LAPACKE_INCLUDE_DIRS ${LAPACKE_INCLUDE_DIR})
ELSE(LAPACKE_FOUND)
  SET(LAPACKE_LIBRARIES)
  SET(LAPACKE_INCLUDE_DIRS)
ENDIF(LAPACKE_FOUND)

MARK_AS_ADVANCED(LAPACKE_INCLUDE_DIR LAPACKE_LIBRARY)
