message(STATUS "CMakeLists.txt in src/")

# Suck up all of the .cc files for the source.  This isn't good CMAKE
# practice, but it makes maintaining this file easier.
file(GLOB dune_tms_source *.cpp)
file(GLOB dune_tms_includes *.h)
# ^^^ Shamelessly ripped from the edep-sim build system, thanks based Clark
#     I'd've written it worse any way

set(source ${dune_tms_source})
set(includes ${dune_tms_includes})

# Build the library.
add_library(TMS_Prod SHARED ${source})

target_include_directories(TMS_Prod PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/toml11/>" # TODO: No chance this is the most correct way
  "$<BUILD_INTERFACE:${ROOT_INCLUDES}>"
  "$<BUILD_INTERFACE:${Geant4_INCLUDES}>"
  "$<BUILD_INTERFACE:${EDEP_INCLUDES}>"
  "$<INSTALL_INTERFACE:include>")

target_link_libraries(TMS_Prod PUBLIC ${EDEP_LIBRARIES} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# Install the TMS_Prod.so library
install(TARGETS TMS_Prod
  EXPORT DUNETMSTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include)

# Install the header files.
install(FILES ${includes} DESTINATION include)
