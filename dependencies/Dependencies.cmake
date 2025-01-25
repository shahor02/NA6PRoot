# Based on O2/dependencies/O2Dependencies.cmake:

# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

include_guard()
include(FeatureSummary)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

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

message(WITH_FLUKA="${WITH_FLUKA}")
set(FlukaRequirement OPTIONAL)
if(DEFINED WITH_FLUKA AND WITH_FLUKA)
  set (FlukaRequirement REQUIRED)  
endif()
find_package(FlukaVMC REQUIRED)
set_package_properties(FlukaVMC PROPERTIES TYPE ${FlukaRequirement})
if(NOT FlukaVMC_FOUND AND "${FlukaRequirement}" STREQUAL "REQUIRED")
  message(FATAL_ERROR "FlukaVMC not found. Please set FlukaVMC_DIR to the directory containing FlukaVMCConfig.cmake.")
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

include_directories(${Geant4_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${VMC_INCLUDE_DIRS})
include_directories(${Geant4VMC_INCLUDE_DIRS})
include_directories($${VGM_INCLUDE_DIRS})
if(FlukaVMC_FOUND)
  include_directories(${FlukaVMC_INCLUDE_DIR})
endif()
