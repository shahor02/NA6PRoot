set(BASE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BASE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Include directoroes
include_directories(${BASE_INCLUDE_DIR})
include_directories(${UTILS_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# Add the library
file(GLOB BASE_SOURCES ${BASE_SRC_DIR}/*.cxx)
file(GLOB BASE_HEADERS ${BASE_INCLUDE_DIR}/*.h)

set(DICTIONARY_HEADERS
  ${BASE_INCLUDE_DIR}/NA6PLayoutParam.h
  ${BASE_INCLUDE_DIR}/MagneticFieldRegion.h
  ${BASE_INCLUDE_DIR}/MagneticField.h
  ) # Add all class headers requiring dictionary
set(DICTIONARY_LINKDEF ${BASE_SRC_DIR}/BaseLinkDef.h)
set(DICTIONARY_OUTPUT BaseDict)

ROOT_GENERATE_DICTIONARY(
    ${DICTIONARY_OUTPUT}
    ${DICTIONARY_HEADERS}
    LINKDEF ${DICTIONARY_LINKDEF}
)

add_library(baseLib SHARED ${BASE_SOURCES} ${BASE_HEADERS} ${DICTIONARY_OUTPUT})

# Set include directories for consumers
target_include_directories(baseLib INTERFACE)

# !!!!!!!!!!!!!!!!!!!!!!
#target_compile_options(baseLib PRIVATE  -ggdb -g -O0)

# Link against required libraries
target_link_libraries(baseLib
    ${ROOT_LIBRARIES}
    fmt::fmt
    FairLogger::FairLogger
    utilsLib
)

# Installation instructions
install(TARGETS baseLib LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES ${BASE_HEADERS} DESTINATION include)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${DICTIONARY_OUTPUT}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${DICTIONARY_OUTPUT}.rootmap DESTINATION lib)
