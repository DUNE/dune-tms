# Hi, this is my first time setting up CMake to build something. Don't judge. :)

cmake_minimum_required(VERSION 3.0)

# Set the project name and version.  The project name will be used to
# import this into other CMAKE files.
project(dune-tms VERSION 0.2.2)

message("-- DUNE TMS Detector Sim and Reconstruction -- ${VERSION}")

# If you got here without git good job honestly
find_package(Git)
#if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.gitmodules")
if(GIT_FOUND AND NOT EXISTS "${PROJECT_SOURCE_DIR}/toml11")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
#        execute_process(COMMAND ls
#                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
#                      message(${GIT_SUBMOD_RESULT})
#execute_process(COMMAND ${GIT_EXECUTALE} submodule update --init  --recursive
        execute_process(COMMAND git submodule update --init  --recursive
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
          message(FATAL_ERROR "${GIT_EXECUTABLE} submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/toml11/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

# Define the options that can be set in the cache, or on the cmake
# command line.
set(CMAKE_BUILD_TYPE Debug)
set(TMS_USE_NEST TRUE CACHE BOOL
  "If true, then make the NEST model available to be directly used")

# Check to see if this is MACOS
if(APPLE)
  message("Apple user: Good luck! I wrote this for Linux so let's hope it works.")
  # If it doesn't just work @liam on slack
  set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif (APPLE)

# Add a target to generate API documentation with Doxygen.  The build
# works fine if this is not found.
find_package(Doxygen)
#if(DOXYGEN_FOUND)
#  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
#    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#  add_custom_target(doc
#    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#    COMMENT "Generating API documentation with Doxygen" VERBATIM)
#else(DOXYGEN_FOUND)
#  add_custom_target(doc
#    COMMENT "Not generating API documentation with Doxygen" VERBATIM)
#endif(DOXYGEN_FOUND)

# Make sure that ROOT is available.  ROOT is absolutely required.
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Look for GEANT4.  If it's found, then build the simulation.
find_package(Geant4)
if (Geant4_FOUND)
  include(${Geant4_USE_FILE})
else (Geant4_FOUND)
  message("GEANT not found")
endif (Geant4_FOUND)

# Check if edep-sim is available
find_package(EDepSim REQUIRED)
if(EDEP_FOUND)
  include(${EDEP_USE_FILE})
endif(EDEP_FOUND)

# Check if toml11 is available, which it should be because submodule
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/toml11/cmake")
#set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/toml11/cmake")
find_package(toml11 REQUIRED)
if(TOML_FOUND)
  include(${TOML_USE_FILE})
endif(TOML_FOUND)


add_subdirectory(src)
add_subdirectory(app)

# Install some extra files.  These are in subdirectories, but are not
# really part of the package.  They are still useful for users to have
# available as examples.
#install(DIRECTORY inputs DESTINATION share/dune-tms
#  FILES_MATCHING
#  PATTERN README*)

#############################################################
#
# Prepare the package so that it can be used with the find_package interface.
#
#############################################################

# Include module with function 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Build the targets description so that the package can be configured
# using find_package.
install(EXPORT DUNETMSTargets
  NAMESPACE TMS::
  DESTINATION lib/cmake/dune-tms)

# Write the 'DUNETMSConfigVersion.cmake' file which can be used to
# check if a version meets the requested properties.
write_basic_package_version_file(
  DUNETMSConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)

# Write the 'DUNETMSConfig.cmake' file so that a user package can
# access this with find_package.
configure_package_config_file(
  PackageConfig.cmake.in
  DUNETMSConfig.cmake
  PATH_VARS CMAKE_INSTALL_PREFIX
  INSTALL_DESTINATION lib/cmake/dune-tms)

# Install the config files.
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/DUNETMSConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/DUNETMSConfigVersion.cmake
  DESTINATION lib/cmake/dune-tms)
