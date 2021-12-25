#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

#-----------------------------------------------------------------------
# - Project definition and basic configuration
#   Version handled manually as project(... VERSION ...) is not used
#   in tests/examples which are subprojects. All calls must use
#   the same form.
project(StandaloneGBP)
set(${PROJECT_NAME}_VERSION_MAJOR  2)
set(${PROJECT_NAME}_VERSION_MINOR  2)
set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build Standalone MC for GBP with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(StandaloneGBP ${PROJECT_SOURCE_DIR}/src/StandaloneGBP.cc ${sources} ${headers})
target_link_libraries(StandaloneGBP ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Configure the installation prefix to allow both local as system-wide installation
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX
        "${PROJECT_SOURCE_DIR}"
        CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build StandaloneGBP. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
include_directories(${PROJECT_SOURCE_DIR}/macros)
set(SCRIPTS
  allpix_check.mac
  gammaconversion.mac
  gaussianBeam.mac
  importBeam.mac
  share.mac
  vis.mac
  )

foreach(_script ${SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/macros/${_script}
    ${PROJECT_BINARY_DIR}/macros/${_script}
    COPYONLY
    )
endforeach()

# add_custom_command(TARGET StandaloneGBP POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy_if_different
#   ${PROJECT_SOURCE_DIR}/macros/*.mac
#   $<TARGET_FILE_DIR:StandaloneGBP>
# )

#----------------------------------------------------------------------------
# Install the executable to 'install' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS StandaloneGBP DESTINATION install)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/macros/ DESTINATION install FILES_MATCHING PATTERN "*.mac")

