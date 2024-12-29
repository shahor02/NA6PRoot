#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(NA6PROOT)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define default installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Default installation directory" FORCE)
endif()

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

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(fmt REQUIRED)
find_package(FairLogger REQUIRED)
find_package(VMC REQUIRED)
find_package(Geant4VMC REQUIRED)
find_package(VGM REQUIRED)
find_package(HepMC3 REQUIRED)

# Locate sources and headers for this project
#
# Define include and source directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${Geant4_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${VMC_INCLUDE_DIRS})
include_directories(${Geant4VMC_INCLUDE_DIRS})
include_directories($${VGM_INCLUDE_DIRS})

# Add subdirectories for each module and prepare paths to its headers
add_subdirectory(utils) # must stay before others
set(UTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils/include)

add_subdirectory(base)
set(BASE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/base/include)

add_subdirectory(sim)
set(SIM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sim/include)

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries
#
add_executable(na6psim NA6PSim.cxx ${TOOLS_FORTRAN_OBJECTS})
target_link_libraries(na6psim
  utilsLib
  baseLib
  simLib
  ${Geant4_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${ROOT_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${VMC_LIBRARIES}
  ${Geant4VMC_LIBRARIES}
  VGM::BaseVGM
  Boost::boost
  Boost::program_options
  fmt::fmt
  FairLogger::FairLogger
  HepMC3::HepMC3
  )
# Explicit dependency on module headers
target_include_directories(na6psim PRIVATE
  ${PROJECT_SOURCE_DIR}/utils/include
  ${PROJECT_SOURCE_DIR}/base/include
  ${PROJECT_SOURCE_DIR}/sim/include
)

# Install instructions (optional)
install(TARGETS na6psim RUNTIME DESTINATION bin)

set(DATA_DIR ${CMAKE_SOURCE_DIR}/data)
# Install all files in the data directory
install(DIRECTORY ${DATA_DIR}/
        DESTINATION share/data
        FILES_MATCHING PATTERN "*.*")

set(TESTS_DIR ${CMAKE_SOURCE_DIR}/test)
# Install all files in the test directory
install(DIRECTORY ${TESTS_DIR}/
        DESTINATION share/test
        FILES_MATCHING PATTERN "*.*")


# Define the content of the init.sh script
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/init.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/init.sh
    @ONLY
)

# Install the generated script to the installation directory
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/init.sh
  DESTINATION .
  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
              GROUP_EXECUTE GROUP_READ
              WORLD_EXECUTE WORLD_READ)

target_compile_options(na6psim PUBLIC -ggdb -g -O0)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ggdb -g")

