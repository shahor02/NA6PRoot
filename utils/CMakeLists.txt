set(UTILS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(UTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Include directoroes
include_directories(${UTILS_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# Add the library
file(GLOB UTILS_SOURCES ${UTILS_SRC_DIR}/*.cxx)
file(GLOB UTILS_HEADERS ${UTILS_INCLUDE_DIR}/*.h)

set(DICTIONARY_HEADERS
  ${UTILS_INCLUDE_DIR}/StringUtils.h
  ${UTILS_INCLUDE_DIR}/ConfigurableParamHelper.h
  ${UTILS_INCLUDE_DIR}/KeyValParam.h
  ${UTILS_INCLUDE_DIR}/MiscUtils.h
  ) # Add all class headers requiring dictionary
set(DICTIONARY_LINKDEF ${UTILS_SRC_DIR}/UtilsLinkDef.h)
set(DICTIONARY_OUTPUT UtilsDict)
set(DICTIONARY_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/dict)

ROOT_GENERATE_DICTIONARY(
    ${DICTIONARY_OUTPUT}
    ${DICTIONARY_HEADERS}
    LINKDEF ${DICTIONARY_LINKDEF}
)

add_library(utilsLib SHARED ${UTILS_SOURCES} ${UTILS_HEADERS} ${DICTIONARY_OUTPUT})

# Set include directories for consumers
target_include_directories(utilsLib INTERFACE)

# Link against required libraries
target_link_libraries(utilsLib
    ${ROOT_LIBRARIES}
    Boost::program_options
    fmt::fmt
    FairLogger::FairLogger
)

# Installation instructions
install(TARGETS utilsLib LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES ${UTILS_HEADERS} DESTINATION include)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${DICTIONARY_OUTPUT}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${DICTIONARY_OUTPUT}.rootmap DESTINATION lib)
