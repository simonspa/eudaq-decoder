########################################################
# CMake file for the pXar framework
CMAKE_MINIMUM_REQUIRED(VERSION 2.4 FATAL_ERROR)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)
########################################################

# Project name
PROJECT( eudecoder )

# Setup the build environment
LIST ( APPEND CMAKE_CXX_FLAGS "-fPIC -O2 -std=c++11" )
LIST ( APPEND CMAKE_LD_FLAGS "-fPIC -O2" )

# Additional packages to be searched for by cmake
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

SET(INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories")
SET(CMAKE_INSTALL_PREFIX "${INSTALL_PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

SET(SOURCE_FILES
  "src/eudaqevd4.cc"
  # "src/main.cc"
  )

# this UI and its libraries depend on ROOT:
FIND_PACKAGE(EUDAQ REQUIRED)
FIND_PACKAGE(ROOT REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM include ${EUDAQ_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${EUDAQ_LIBRARY} ${ROOT_LIBRARIES})
