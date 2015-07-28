# - Try to find pxarCore used in native reader library for exporting RAW to pxarCore format needed for reconstruction
# Once done this will define
#  EUDAQ_FOUND - System has EUDAQ
#  EUDAQ_INCLUDE_DIR - The EUDAQ main include directories
#  EUDAQ_LIBRARY - The libraries needed to use EUDAQ

MESSAGE(STATUS "Looking for EUDAQ...")

find_path(EUDAQ_INCLUDE_DIR eudaq/FileReader.hh
  HINTS "${EUDAQPATH}/main/include" "$ENV{EUDAQPATH}/main/include")

find_library(EUDAQ_LIBRARY NAMES EUDAQ
  HINTS "${EUDAQPATH}/lib" "$ENV{EUDAQPATH}/lib")

MESSAGE(STATUS "include: ${EUDAQ_INCLUDE_DIR} ${EUDAQ_LIBRARY}")

#include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EUDAQ_FOUND to TRUE
# if all listed variables are TRUE
#find_package_handle_standard_args(pxarCore
#  REQUIRED_VARS EUDAQ_LIBRARY EUDAQ_API_INCLUDE_DIR EUDAQ_UTILS_INCLUDE_DIR PXAR_UTIL_INCLUDE_DIR)

IF(EUDAQ_LIBRARY AND EUDAQ_INCLUDE_DIR)
   SET(EUDAQ_FOUND TRUE)
   MESSAGE(STATUS "Found EUDAQ library and headers.")
   SET(EUDAQ_INCLUDE_DIR "$ENV{EUDAQPATH}/main/include/")
   MESSAGE(STATUS "include: ${EUDAQ_INCLUDE_DIR} ${EUDAQ_LIBRARY}")
ENDIF(EUDAQ_LIBRARY AND EUDAQ_INCLUDE_DIR)

mark_as_advanced(EUDAQ_LIBRARY EUDAQ_INCLUDE_DIR)