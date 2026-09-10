message(STATUS "CMakeLists.txt in src/")

# Suck up all of the .cc files for the source, including subdirectories like
# Cluster3D/ (space-point building, clustering, and pattern-recognition code).
# This isn't good CMAKE practice, but it makes maintaining this file easier.
file(GLOB_RECURSE dune_tms_source *.cpp)
file(GLOB_RECURSE dune_tms_includes *.h)
# ^^^ Shamelessly ripped from the edep-sim build system, thanks based Clark
#     I'd've written it worse any way

# attic/ is a graveyard for retired code -- not meant to build (its own
# include is already broken: TMS_Classes.cpp asks for "TMSClasses.h", which
# doesn't exist). Keep it out of the recursive glob above.
list(FILTER dune_tms_source EXCLUDE REGEX "/attic/")
list(FILTER dune_tms_includes EXCLUDE REGEX "/attic/")

set(source ${dune_tms_source})
set(includes ${dune_tms_includes})

# Build the library.
add_library(TMS_Prod SHARED ${source})

# Enable compiler warnings to catch issues early (GCC/Clang only for portability)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(TMS_Prod PRIVATE -Wall -Wextra)
endif()

target_include_directories(TMS_Prod PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Cluster3D>"
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/toml11/>" # TODO: No chance this is the most correct way
  "$<BUILD_INTERFACE:${ROOT_INCLUDES}>"
  "$<BUILD_INTERFACE:${Geant4_INCLUDES}>"
  "$<BUILD_INTERFACE:${EDEP_INCLUDES}>"
  "$<INSTALL_INTERFACE:include>")

target_link_libraries(TMS_Prod PUBLIC ${EDEP_LIBRARIES} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ROOT::EG)

# Install the TMS_Prod.so library
install(TARGETS TMS_Prod
  EXPORT DUNETMSTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include)

# Install the header files.
install(FILES ${includes} DESTINATION include)
