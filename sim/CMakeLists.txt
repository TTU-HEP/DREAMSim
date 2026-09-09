#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(CaloX)

#----------------------------------------------------------------------------
# Default to an optimised build.  Compiling takes about the same either way (a
# few seconds), while an unoptimised binary runs about twice as slow: 10.6 s per
# 100 GeV e+ event against 5.0 s.  Override with -DCMAKE_BUILD_TYPE=Debug when
# stepping through the code.
#
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "CaloX build type: ${CMAKE_BUILD_TYPE}")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#~~~~~~~~~~~~~~~~~~~~
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS})
#~~~~~~~~~~~~~~~~~~~~


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
add_executable(CaloXSim CaloXSim.cc ${sources} ${headers})
target_link_libraries(CaloXSim ${Geant4_LIBRARIES})
target_link_libraries(CaloXSim ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build CaloX. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLECALOX_SCRIPTS
  vis.mac
  paramBatch03_single.mac
  paramBatch03_single_photon.mac
  paramBatch03_single_rotated.mac
  paramBatch03_jet.mac
  )


foreach(_script ${EXAMPLECALOX_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Copy the fiber map next to the executable, so that the default
# fiberMapFile of data/fibermap.json resolves from the build directory.
#
configure_file(
  ${PROJECT_SOURCE_DIR}/data/fibermap.json
  ${PROJECT_BINARY_DIR}/data/fibermap.json
  COPYONLY
  )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CaloXSim DESTINATION bin)
