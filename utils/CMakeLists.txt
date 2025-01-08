set(THIS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(THIS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(MOD_NAME "utils")
set(LIB_NAME "${MOD_NAME}Lib")

# Include directoroes
include_directories(${THIS_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

file(GLOB THIS_SOURCES ${THIS_SRC_DIR}/*.cxx)
file(GLOB THIS_HEADERS ${THIS_INCLUDE_DIR}/*.h)

# Add the library
add_library(utilsLib SHARED ${THIS_SOURCES} ${THIS_HEADERS})

set(DICTIONARY_HEADERS
  StringUtils.h
  ConfigurableParamHelper.h
  KeyValParam.h
  MiscUtils.h
  ) # Add all class headers requiring dictionary

set(DICTIONARY_LINKDEF ${THIS_SRC_DIR}/${MOD_NAME}LinkDef.h)
set(DICTIONARY_OUTPUT ${LIB_NAME}Dict)

ROOT_GENERATE_DICTIONARY(
    ${DICTIONARY_OUTPUT}
    ${DICTIONARY_HEADERS}
    LINKDEF ${DICTIONARY_LINKDEF}
    MODULE ${LIB_NAME}
)

# Set include directories for consumers
target_include_directories(${LIB_NAME} INTERFACE ${THIS_INCLUDE_DIR})

# !!!!!!!!!!!!!!!!!!!!!!
#target_compile_options(${LIB_NAME} PRIVATE -ggdb -g -O0)

# Link against required libraries
target_link_libraries(${LIB_NAME}
    ${ROOT_LIBRARIES}
    Boost::program_options
    fmt::fmt
    FairLogger::FairLogger
)

# Installation instructions
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES ${THIS_HEADERS} DESTINATION include)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}.rootmap DESTINATION lib)
