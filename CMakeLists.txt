cmake_minimum_required (VERSION 3.5)

project (PATSMA VERSION 1.0.0)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(Green       "${Esc}[32m")
  set(ColourReset "${Esc}[m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldCyan    "${Esc}[1;36m")
endif()

function(message)
  list(GET ARGV 0 MessageType)
  if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
    list(REMOVE_AT ARGV 0)
    _message(${MessageType} "${BoldRed}${ARGV}${ColourReset}")
  elseif(MessageType STREQUAL WARNING)
    list(REMOVE_AT ARGV 0)
    _message(STATUS "${BoldYellow}${ARGV}${ColourReset}")
  elseif(MessageType STREQUAL AUTHOR_WARNING)
    list(REMOVE_AT ARGV 0)
    _message(STATUS "${BoldCyan}${ARGV}${ColourReset}")
  elseif(MessageType STREQUAL STATUS)
    list(REMOVE_AT ARGV 0)
    _message(STATUS "${Green}${ARGV}${ColourReset}")
  else()
    _message("${ARGV}")
  endif()
endfunction()

option (VERBOSE "Enable verbose mode." OFF)
option (TEST "Compile tests." OFF)
option (OPENMP "Enable OpenMP pragma." ON)

if (OPENMP)
  find_package (OpenMP QUIET)
  if (NOT OPENMP_FOUND)
    message (FATAL_ERROR "OpenMP must exist to compile FWI when OpenMP is enabled.")
  else ()
    message (STATUS "OpenMP ${OpenMP_CXX_VERSION} found & enabled.")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  endif ()
endif ()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message (STATUS "Setting build type to 'Release' as none was specified.")
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
else ()
  STRING (TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)
  set (DEBUG OFF)
  if (CMAKE_BUILD_TYPE_LOWER STREQUAL "debug")
    set (DEBUG ON)
  endif ()
endif ()


# Enable code testing
if(TEST)
  message (STATUS "Setting build with TESTING.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTEST -Wall")
endif()

# Enable verbosing
if (VERBOSE)
  message (STATUS "Setting build with VERBOSE.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE -Wall")
endif ()

if (DEBUG)
  message (STATUS "Setting build type to 'Debug'.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DDEBUG -Wall")
endif ()

## LIBRARY INSTALL CONFIGURATION ##
# Gather source and header files using 'file' command
file(GLOB TEMPLTS "src/*.tpp")
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "src/*.hpp")

# Create a static library from the source files
add_library(${PROJECT_NAME}_lib STATIC ${SOURCES})

# Specify include directories, compile options, and libraries for the library
target_include_directories(${PROJECT_NAME}_lib PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_compile_options(${PROJECT_NAME}_lib PUBLIC ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(${PROJECT_NAME}_lib ${MPI_CXX_LIBRARIES})

# Install headers and the library
install(FILES ${HEADERS} ${TEMPLTS} DESTINATION include/)
install(TARGETS ${PROJECT_NAME}_lib DESTINATION lib/)

# Create an uppercase version of the project name
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)

# Include CMakePackageConfigHelpers to generate package config files
include(CMakePackageConfigHelpers)

# Generate the version file for the package
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_libConfigVersion.cmake"
  VERSION 1.0
  COMPATIBILITY AnyNewerVersion
)

# Generate the package config file
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_libConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake/
)

# Install the generated package config files
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_libConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_libConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake/
)