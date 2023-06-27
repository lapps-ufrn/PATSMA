cmake_minimum_required (VERSION 3.5)

project (autotuning)

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

separate_arguments(MPI_CXX_COMPILE_FLAGS)

# Enable code testing
if(TEST)
  message (STATUS "Setting build with TESTING.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTEST")
endif()

# Enable verbosing
if (VERBOSE)
  message (STATUS "Setting build with VERBOSE.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
endif ()

if (DEBUG)
  message (STATUS "Setting build type to 'Debug'.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DDEBUG")
endif ()

# # Find source files
# file(GLOB SOURCES src/*.cpp)
# message(${SOURCES})

# # Create shared library
# add_library(${PROJECT_NAME} SHARED ${SOURCES})

# Install library
# install(TARGETS ${PROJECT_NAME} DESTINATION examples/AT_LIB/${PROJECT_NAME}/lib)

# # Install library headers
# file(GLOB HEADERS src/*.h)
# target_include_directories (${PROJECT_NAME} PUBLIC ${MPI_CXX_INCLUDE_DIRS})
# target_compile_options (${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS})

# install(FILES ${HEADERS} DESTINATION examples/AT_LIB/${PROJECT_NAME}/include)
# install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Config
#         ARCHIVE  DESTINATION examples/AT_LIB/${PROJECT_NAME}
#         LIBRARY  DESTINATION examples/AT_LIB/${PROJECT_NAME}
#         RUNTIME  DESTINATION examples/AT_LIB/${PROJECT_NAME})

# install(EXPORT ${PROJECT_NAME}Config DESTINATION examples/AT_LIB/${PROJECT_NAME})

# export(TARGETS ${PROJECT_NAME} FILE examples/AT_LIB/${PROJECT_NAME}/${PROJECT_NAME}Config.cmake)

add_subdirectory (src)

# Project test source files
if (TEST)
  add_subdirectory (test)
endif ()