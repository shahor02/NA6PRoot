set(SIM_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SIM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Include directoroes
include_directories(${SIM_INCLUDE_DIR})
include_directories(${UTILS_INCLUDE_DIR})
include_directories(${BASE_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${VMC_INCLUDE_DIRS})

file(GLOB SIM_SOURCES ${SIM_SRC_DIR}/*.cxx)
file(GLOB SIM_HEADERS ${SIM_INCLUDE_DIR}/*.h)

set(DICTIONARY_HEADERS
  NA6PMCStack.h
  NA6PGenerator.h
  NA6PGenCutParam.h
  NA6PMCEventHeader.h
  NA6PMCGenHeader.h
  NA6PGenBox.h
  NA6PBaseHit.h
  NA6PVerTelHit.h
  NA6PMuonSpecHit.h
  ) # Add all class headers requiring dictionary
set(DICTIONARY_LINKDEF ${SIM_SRC_DIR}/SimLinkDef.h)
set(DICTIONARY_OUTPUT SimDict)

ROOT_GENERATE_DICTIONARY(
    ${DICTIONARY_OUTPUT}
    ${DICTIONARY_HEADERS}
    LINKDEF ${DICTIONARY_LINKDEF}
)

# Add the library
add_library(simLib SHARED ${SIM_SOURCES} ${SIM_HEADERS} ${DICTIONARY_OUTPUT})

# Set include directories for consumers
target_include_directories(simLib INTERFACE ${SIM_INCLUDE_DIR})

# !!!!!!!!!!!!!!!!!!!!!!
#target_compile_options(simLib PRIVATE -ggdb -g -O0)


# Link against required libraries
target_link_libraries(simLib
    ${ROOT_LIBRARIES}
    fmt::fmt
    FairLogger::FairLogger
    ${VMC_LIBRARIES}
    baseLib
    utilsLib
)

# Installation instructions
install(TARGETS simLib LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES ${SIM_HEADERS} DESTINATION include)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${DICTIONARY_OUTPUT}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${DICTIONARY_OUTPUT}.rootmap DESTINATION lib)
