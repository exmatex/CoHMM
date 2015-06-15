# - Find cblas
# Find the native CBLAS headers and libraries.
#
#  CBLAS_INCLUDE_DIRS - where to find cblas.h, etc.
#  CBLAS_LIBRARIES    - List of libraries when using cblas.
#  CBLAS_FOUND        - True if cblas found.

# Look for the header file.
FIND_PATH(CBLAS_INCLUDE_DIR NAMES cblas.h PATH_SUFFIXES cblas)

# Look for the library.
FIND_LIBRARY(CBLAS_LIBRARY NAMES cblas refcblas)

#cblas needs blas
FIND_PACKAGE(BLAS)

# handle the QUIETLY and REQUIRED arguments and set CBLAS_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CBLAS DEFAULT_MSG CBLAS_LIBRARY BLAS_LIBRARIES CBLAS_INCLUDE_DIR)

# Copy the results to the output variables.
IF(CBLAS_FOUND)
  SET(CBLAS_LIBRARIES ${CBLAS_LIBRARY} ${BLAS_LIBRARIES})
  SET(CBLAS_INCLUDE_DIRS ${CBLAS_INCLUDE_DIR})
ELSE(CBLAS_FOUND)
  SET(CBLAS_LIBRARIES)
  SET(CBLAS_INCLUDE_DIRS)
ENDIF(CBLAS_FOUND)

MARK_AS_ADVANCED(CBLAS_INCLUDE_DIR CBLAS_LIBRARY)
